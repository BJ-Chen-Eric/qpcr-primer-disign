#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import re
import os
import subprocess
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

import requests

ENSEMBL = "https://rest.ensembl.org"
SPECIES = "homo_sapiens"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


@dataclass
class Exon:
    exon_id: str
    start: int
    end: int
    strand: int


@dataclass
class JunctionDesign:
    gene: str
    transcript_id: str
    transcript_version: Optional[int]
    exon_left: str
    exon_right: str
    exon_left_num: int
    exon_right_num: int
    left_primer: str
    right_primer: str
    left_tm: float
    right_tm: float
    pair_tm: Optional[float]
    pair_penalty: Optional[float]
    product_size: int
    template: str
    junction_pos: int
    left_pos: int
    left_len: int
    right_pos: int
    right_len: int
    which_spans: str  # LEFT / RIGHT
    junction_overlap_left: int
    junction_overlap_right: int
    score: float


def ensembl_get(path: str, params=None):
    url = f"{ENSEMBL}{path}"
    r = requests.get(url, headers=HEADERS, params=params, timeout=60)
    if not r.ok:
        raise RuntimeError(f"Ensembl request failed {r.status_code}: {r.text[:2000]}")
    return r.json()


def fetch_gene_with_transcripts(gene_symbol: str, species: str):
    """
    Try symbol lookup first; if not found, fall back to ID lookup.
    """
    try:
        return ensembl_get(f"/lookup/symbol/{species}/{gene_symbol}", params={"expand": "1"})
    except RuntimeError:
        return ensembl_get(f"/lookup/id/{gene_symbol}", params={"expand": "1"})


def pick_transcript(gene_json: dict) -> dict:
    txs = gene_json.get("Transcript", [])
    if not txs:
        raise RuntimeError("找不到 transcripts（gene symbol 可能不對或 Ensembl 無資料）。")

    canon = [t for t in txs if t.get("is_canonical") == 1]
    if canon:
        return canon[0]

    def tx_score(t):
        ex = t.get("Exon", [])
        return sum(abs(e["end"] - e["start"]) + 1 for e in ex)

    return sorted(txs, key=tx_score, reverse=True)[0]


def exons_in_tx(tx: dict) -> List[Exon]:
    ex = tx.get("Exon", [])
    if not ex:
        raise RuntimeError("此 transcript 沒有 exon 資訊。")

    strand = int(tx["strand"])
    ex_sorted = sorted(ex, key=lambda e: e["start"])
    if strand == -1:
        ex_sorted = list(reversed(ex_sorted))
    return [Exon(e["id"], int(e["start"]), int(e["end"]), strand) for e in ex_sorted]


def fetch_exon_seq(exon_id: str) -> str:
    j = ensembl_get(f"/sequence/id/{exon_id}")
    seq = j.get("seq")
    if not seq:
        raise RuntimeError(f"exon {exon_id} 沒有回傳序列。")
    return seq.upper()


def fetch_transcript_cdna(transcript_id: str) -> str:
    j = ensembl_get(f"/sequence/id/{transcript_id}", params={"type": "cdna"})
    seq = j.get("seq")
    if not seq:
        raise RuntimeError(f"transcript {transcript_id} 沒有回傳 cDNA 序列。")
    return seq.upper()


def build_junction_template(exon_left_seq: str, exon_right_seq: str, flank: int) -> Tuple[str, int]:
    left = exon_left_seq[-flank:] if len(exon_left_seq) > flank else exon_left_seq
    right = exon_right_seq[:flank] if len(exon_right_seq) > flank else exon_right_seq
    template = (left + right).upper()
    junction_pos = len(left)
    return template, junction_pos


def parse_region_arg(region: str, tx_len: int) -> Tuple[int, int]:
    """
    region: "start-end" (1-based, inclusive). Return (start0, end0_exclusive).
    """
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", region)
    if not m:
        raise ValueError("region 格式錯誤，請用 start-end（1-based）")
    start = int(m.group(1))
    end = int(m.group(2))
    if start < 1 or end < 1 or start > end:
        raise ValueError("region 需符合 1-based 且 start <= end")
    if end > tx_len:
        raise ValueError(f"region 超出 transcript 長度（{tx_len}）")
    return start - 1, end


def build_region_template(cdna: str, start0: int, end0: int) -> Tuple[str, int]:
    """
    Slice transcript cDNA by 0-based [start0, end0) range.
    Returns (template, template_tx_start).
    """
    return cdna[start0:end0], start0


def build_center_template(cdna: str, flank: int) -> Tuple[str, int]:
    """
    Build a ~2*flank template centered on transcript cDNA.
    Returns (template, template_tx_start).
    """
    n = len(cdna)
    if n <= 2 * flank:
        return cdna, 0
    center = n // 2
    start0 = max(0, center - flank)
    end0 = min(n, start0 + 2 * flank)
    start0 = max(0, end0 - 2 * flank)
    return cdna[start0:end0], start0


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    gc = sum(1 for b in s if b in ("G", "C"))
    return gc * 100.0 / len(s)


def map_template_pos_to_tx(
    pos: int,
    left_exon_len: int,
    right_exon_len: int,
    left_exon_tx_start: int,
    right_exon_tx_start: int,
    left_part_len: int,
) -> int:
    """
    將 template 上的 0-based 位置映射到 transcript (cDNA) 0-based 位置。
    template = left_exon_tail + right_exon_head
    """
    if pos < left_part_len:
        offset_in_left = (left_exon_len - left_part_len) + pos
        return left_exon_tx_start + offset_in_left
    return right_exon_tx_start + (pos - left_part_len)


def map_template_pos_to_tx_offset(pos: int, template_tx_start: int) -> int:
    return template_tx_start + pos


def choose_junction_indices(n_exons: int, max_shift: int = 4) -> List[int]:
    """
    以「中間 junction」為中心，左右擴張嘗試多個 junction。
    回傳的是 left exon 的 index（junction 在 i -> i+1）。
    """
    if n_exons < 2:
        return []
    center = max(0, n_exons // 2 - 1)
    order = [center]
    for s in range(1, max_shift + 1):
        if center - s >= 0:
            order.append(center - s)
        if center + s < n_exons - 1:
            order.append(center + s)
    return order


def run_primer3_pairs(
    primer3_core: str,
    template: str,
    prod_min: int,
    prod_max: int,
    num_return: int,
) -> Dict[str, str]:
    """
    讓 primer3 產生多組 primer pairs，回傳 key=value 字典。
    """
    p3 = f"""\
SEQUENCE_ID=JUNC
SEQUENCE_TEMPLATE={template}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=57.0
PRIMER_MAX_TM=63.0
PRIMER_PRODUCT_SIZE_RANGE={prod_min}-{prod_max}
PRIMER_NUM_RETURN={num_return}
PRIMER_MAX_POLY_X=5
PRIMER_GC_CLAMP=1
=
"""
    proc = subprocess.run(
        [primer3_core],
        input=p3.encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.decode(errors="replace")[:2000])

    kv: Dict[str, str] = {}
    for line in proc.stdout.decode(errors="replace").splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()
    return kv


def parse_pair(
    kv: Dict[str, str], i: int
) -> Optional[Tuple[str, int, int, float, str, int, int, float, Optional[float], Optional[float], int]]:
    """
    解析第 i 組 primer pair。
    回傳 (Lseq, Lpos, Llen, Ltm, Rseq, Rpos, Rlen, Rtm, pair_tm, pair_penalty, product_size)
    其中 PRIMER_LEFT_i = "pos,len"
         PRIMER_RIGHT_i = "pos,len"  (pos 是右引子的 5' 位置在 template 上)
    """
    lseq = kv.get(f"PRIMER_LEFT_{i}_SEQUENCE")
    rseq = kv.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
    lpl = kv.get(f"PRIMER_LEFT_{i}")
    rpl = kv.get(f"PRIMER_RIGHT_{i}")
    ltm = kv.get(f"PRIMER_LEFT_{i}_TM")
    rtm = kv.get(f"PRIMER_RIGHT_{i}_TM")
    ptm = kv.get(f"PRIMER_PAIR_{i}_PRODUCT_TM")
    ppen = kv.get(f"PRIMER_PAIR_{i}_PENALTY")
    psz = kv.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE")

    if not (lseq and rseq and lpl and rpl and ltm and rtm and psz):
        return None

    lpos, llen = [int(x) for x in lpl.split(",")]
    rpos, rlen = [int(x) for x in rpl.split(",")]
    pair_tm = float(ptm) if ptm else None
    pair_penalty = float(ppen) if ppen else None
    return lseq, lpos, llen, float(ltm), rseq, rpos, rlen, float(rtm), pair_tm, pair_penalty, int(psz)


def pick_metric(kv: Dict[str, str], keys: List[str]) -> Optional[str]:
    for k in keys:
        v = kv.get(k)
        if v is not None and v != "":
            return v
    return None

def get_float_metric(kv: Dict[str, str], keys: List[str]) -> Optional[float]:
    v = pick_metric(kv, keys)
    if v is None:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def spans_junction(start: int, length: int, junction_pos: int) -> bool:
    end = start + length
    return start < junction_pos < end

def junction_overlap_bp(start: int, length: int, junction_pos: int) -> int:
    """
    計算 primer 與 junction 的重疊 bp 數。
    primer 覆蓋區間為 [start, start+length)；junction_pos 是切點（在兩段之間）。
    overlap = 左側覆蓋 bp + 右側覆蓋 bp
    若沒有跨 junction，回傳 0。
    """
    end = start + length
    if not (start < junction_pos < end):
        return 0
    left_overlap = junction_pos - start
    right_overlap = end - junction_pos
    return min(left_overlap, length) + min(right_overlap, length)


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))

def gc_penalty(gc: float, lo: float = 35.0, hi: float = 65.0) -> float:
    """
    GC 在 [lo, hi] 內 -> 0 penalty
    超出 -> 線性懲罰（可自行調）
    """
    if gc < lo:
        return (lo - gc) / lo
    if gc > hi:
        return (gc - hi) / (100.0 - hi)
    return 0.0

def size_penalty(product_size: int, opt: int) -> float:
    """
    產物長度偏離 opt 的懲罰（相對偏差）
    """
    if opt <= 0:
        return 0.0
    return abs(product_size - opt) / opt

def score_pair(
    kv: Dict[str, str],
    i: int,
    lseq: str,
    lpos: int,
    llen: int,
    ltm: float,
    rseq: str,
    rpos: int,
    rlen: int,
    rtm: float,
    pair_tm: Optional[float],
    product_size: int,
    junction_pos: int,
    opt_tm: float = 60.0,
    opt_size: int = 110,
) -> Tuple[float, Dict[str, float]]:
    """
    回傳 (score, breakdown)。
    breakdown 方便 debug/調參。
    """
    # GC
    lgc = gc_percent(lseq)
    rgc = gc_percent(rseq)

    # primer3 metrics（取 TH 優先，沒有就用非 TH）
    l_self_any = get_float_metric(kv, [f"PRIMER_LEFT_{i}_SELF_ANY_TH", f"PRIMER_LEFT_{i}_SELF_ANY"]) or 0.0
    l_self_end = get_float_metric(kv, [f"PRIMER_LEFT_{i}_SELF_END_TH", f"PRIMER_LEFT_{i}_SELF_END"]) or 0.0
    l_hairpin  = get_float_metric(kv, [f"PRIMER_LEFT_{i}_HAIRPIN_TH", f"PRIMER_LEFT_{i}_HAIRPIN"]) or 0.0

    r_self_any = get_float_metric(kv, [f"PRIMER_RIGHT_{i}_SELF_ANY_TH", f"PRIMER_RIGHT_{i}_SELF_ANY"]) or 0.0
    r_self_end = get_float_metric(kv, [f"PRIMER_RIGHT_{i}_SELF_END_TH", f"PRIMER_RIGHT_{i}_SELF_END"]) or 0.0
    r_hairpin  = get_float_metric(kv, [f"PRIMER_RIGHT_{i}_HAIRPIN_TH", f"PRIMER_RIGHT_{i}_HAIRPIN"]) or 0.0

    p_compl_any = get_float_metric(kv, [f"PRIMER_PAIR_{i}_COMPL_ANY_TH", f"PRIMER_PAIR_{i}_COMPL_ANY"]) or 0.0
    p_compl_end = get_float_metric(kv, [f"PRIMER_PAIR_{i}_COMPL_END_TH", f"PRIMER_PAIR_{i}_COMPL_END"]) or 0.0

    # penalties
    tm_off_L = abs(ltm - opt_tm)
    tm_off_R = abs(rtm - opt_tm)
    tm_diff  = abs(ltm - rtm)

    gc_off = gc_penalty(lgc) + gc_penalty(rgc)
    sz_off = size_penalty(product_size, opt_size)

    # 權重（你可依 lab 偏好調整）
    w_tm_off = 0.8
    w_tm_diff = 1.2
    w_gc = 2.0
    w_size = 0.7

    w_self_any = 0.25
    w_self_end = 0.60
    w_hairpin  = 0.08

    w_pair_any = 0.35
    w_pair_end = 0.60

    penalty = 0.0
    penalty += w_tm_off * (tm_off_L + tm_off_R)
    penalty += w_tm_diff * tm_diff
    penalty += w_gc * gc_off
    penalty += w_size * (sz_off * 10.0)  # 放大一下 size 的量級（可調）

    penalty += w_self_any * (l_self_any + r_self_any)
    penalty += w_self_end * (l_self_end + r_self_end)
    penalty += w_hairpin  * (l_hairpin  + r_hairpin)

    penalty += w_pair_any * p_compl_any
    penalty += w_pair_end * p_compl_end

    # score：越大越好
    score = -penalty

    breakdown = {
        "tm_off_sum": tm_off_L + tm_off_R,
        "tm_diff": tm_diff,
        "gc_pen": gc_off,
        "size_pen": sz_off,
        "l_self_any": l_self_any,
        "l_self_end": l_self_end,
        "l_hairpin": l_hairpin,
        "r_self_any": r_self_any,
        "r_self_end": r_self_end,
        "r_hairpin": r_hairpin,
        "pair_compl_any": p_compl_any,
        "pair_compl_end": p_compl_end,
        "lgc": lgc,
        "rgc": rgc,
    }
    return score, breakdown


def main():
    parser = argparse.ArgumentParser(description="設計 exon–exon junction 的 SYBR RT-qPCR primers（human / Ensembl）")
    parser.add_argument("--gene", help="Gene symbol（例如 TP53、IL6、CXCL8）")
    parser.add_argument("--species", default=SPECIES, help="Ensembl species（例如 homo_sapiens, mus_musculus）")
    parser.add_argument("--sequence", help="User-provided DNA sequence (A/C/G/T/N), no FASTA header")
    parser.add_argument("--primer3", required=True, help="primer3_core 可執行檔的絕對路徑")
    parser.add_argument("--prod-min", type=int, default=70, help="產物最小長度（bp），預設 70")
    parser.add_argument("--prod-max", type=int, default=180, help="產物最大長度（bp），預設 180")
    parser.add_argument("--flank", type=int, default=150, help="junction template 每側長度（bp），預設 150（總長約 300）")
    parser.add_argument("--region", help="指定 cDNA 區段（1-based, inclusive），格式 start-end，例如 200-400")
    parser.add_argument("--allow-single-exon", action="store_true", help="允許單 exon transcript（無 junction）設計")
    parser.add_argument("--num-return", type=int, default=50, help="primer3 回傳幾組 pair，預設 50")
    parser.add_argument("--max-junction-shift", type=int, default=4, help="以中間 junction 為中心，左右最多換幾個 junction，預設 4")
    parser.add_argument("--opt-amplicon", type=int, default=110, help="偏好的產物長度（bp），用於 ranking（預設 110）")
    parser.add_argument("--top", type=int, default=5, help="全局排序後輸出前幾名候選（預設 5）")
    parser.add_argument("--min-candidates-per-junction", type=int, default=3, help="每個 junction 至少保留幾組候選，少於此數則換下一個 junction")
    parser.add_argument("--out", help="輸出結果到 txt 檔案（例如 /abs/path/result.txt）")
    parser.add_argument(
        "--span",
        choices=["ANY", "LEFT", "RIGHT"],
        default="ANY",
        help="指定哪一條 primer 必須跨 junction：ANY（預設，任一條）/ LEFT / RIGHT",
    )
    parser.add_argument(
        "--min-junction-overlap",
        type=int,
        default=5,
        help="跨 junction 的最小重疊 bp 數（預設 5；設 0 表示不限制）",
    )

    args = parser.parse_args()

    gene = args.gene.strip() if args.gene else None
    species = args.species.strip().lower() if args.species else None
    user_sequence = args.sequence.strip() if args.sequence else None
    primer3_core = args.primer3

    if not os.path.isfile(primer3_core):
        raise FileNotFoundError(f"找不到 primer3_core：{primer3_core}")
    if not os.access(primer3_core, os.X_OK):
        raise PermissionError(f"primer3_core 無可執行權限：{primer3_core}")

    if user_sequence:
        seq = re.sub(r"\\s+", "", user_sequence).upper()
        if not re.fullmatch(r"[ACGTN]+", seq):
            raise ValueError("sequence 只允許 A/C/G/T/N")
        gene_name = None
        tx_id = "user_sequence"
        tx_version = None
        cdna_seq = seq
        tx_len = len(seq)
        exons = [Exon("EXON1", 1, tx_len, 1)]
        exon_lengths = [tx_len]
        exon_tx_starts = [0]
    else:
        if not gene:
            raise ValueError("請提供 --gene 或 --sequence")
        gene_json = fetch_gene_with_transcripts(gene, species)
        gene_name = gene_json.get("display_name") or gene_json.get("name")
        tx = pick_transcript(gene_json)
        tx_id = tx["id"]
        tx_version = tx.get("version")
        exons = exons_in_tx(tx)
        exon_lengths = [abs(e.end - e.start) + 1 for e in exons]
        exon_tx_starts: List[int] = []
        acc = 0
        for L in exon_lengths:
            exon_tx_starts.append(acc)
            acc += L
        tx_len = acc

    junction_mode = len(exons) >= 2
    if not junction_mode and not args.allow_single_exon:
        print("Warning: single-exon transcript detected; switching to single-exon mode.", file=os.sys.stderr)
        print(
            "Warning: single-exon assay: cannot avoid gDNA by junction. Do DNase and no-RT control; "
            "in-silico specificity is strongly recommended.",
            file=os.sys.stderr,
        )

    if junction_mode:
        junction_indices = choose_junction_indices(len(exons), max_shift=args.max_junction_shift)
        if not junction_indices:
            raise RuntimeError("此 transcript exon 少於 2 個：無法做 exon–exon junction primer。")
    else:
        junction_indices = []

    # cDNA sequence is only needed for region/single-exon templates
    if not user_sequence:
        cdna_seq = None
        if args.region or not junction_mode:
            cdna_seq = fetch_transcript_cdna(tx_id)

    region_start0 = None
    region_end0 = None
    template_tx_start = None
    region_template = None
    if args.region:
        region_start0, region_end0 = parse_region_arg(args.region, tx_len)
        region_template, template_tx_start = build_region_template(cdna_seq, region_start0, region_end0)

    last_err: Optional[str] = None
    all_candidates = []

    for idx in junction_indices:
        eL = exons[idx]
        eR = exons[idx + 1]
        eL_num = idx + 1
        eR_num = idx + 2
        total_junctions = len(exons) - 1
        eL_len = exon_lengths[idx]
        eR_len = exon_lengths[idx + 1]
        eL_tx_start = exon_tx_starts[idx]
        eR_tx_start = exon_tx_starts[idx + 1]
        junction_tx_pos = eL_tx_start + eL_len

        if args.region:
            template = region_template
            jpos = junction_tx_pos - region_start0
            left_part_len = None
            map_mode = "offset"
            template_start_tx = template_tx_start
        else:
            seqL = fetch_exon_seq(eL.exon_id)
            seqR = fetch_exon_seq(eR.exon_id)
            template, jpos = build_junction_template(seqL, seqR, flank=args.flank)
            left_part_len = jpos
            map_mode = "junction"
            template_start_tx = None

        kv = run_primer3_pairs(
            primer3_core=primer3_core,
            template=template,
            prod_min=args.prod_min,
            prod_max=args.prod_max,
            num_return=args.num_return,
        )
        junction_candidates = []

        # 在所有 pair 中找符合條件者
        for i in range(args.num_return):
            parsed = parse_pair(kv, i)
            if not parsed:
                continue
            lseq, lpos, llen, ltm, rseq, rpos, rlen, rtm, ptm, ppen, psz = parsed

            l_span = spans_junction(lpos, llen, jpos)
            r_span = spans_junction(rpos, rlen, jpos)

            l_ov = junction_overlap_bp(lpos, llen, jpos) if l_span else 0
            r_ov = junction_overlap_bp(rpos, rlen, jpos) if r_span else 0

            # span 規則
            want = args.span
            if want == "LEFT" and not l_span:
                continue
            if want == "RIGHT" and not r_span:
                continue
            if want == "ANY" and not (l_span or r_span):
                continue

            # overlap 門檻
            min_ov = max(0, int(args.min_junction_overlap))
            if min_ov > 0:
                if want == "LEFT":
                    if l_ov < min_ov:
                        continue
                elif want == "RIGHT":
                    if r_ov < min_ov:
                        continue
                else:  # ANY
                    if max(l_ov, r_ov) < min_ov:
                        continue

            # which_spans（若兩條都跨，選 overlap 較大者）
            if l_span and r_span:
                which = "LEFT" if l_ov >= r_ov else "RIGHT"
            else:
                which = "LEFT" if l_span else "RIGHT"

            score, breakdown = score_pair(
                kv=kv, i=i,
                lseq=lseq, lpos=lpos, llen=llen, ltm=ltm,
                rseq=rseq, rpos=rpos, rlen=rlen, rtm=rtm,
                pair_tm=ptm, product_size=psz,
                junction_pos=jpos,
                opt_tm=60.0,
                opt_size=args.opt_amplicon,
            )

            junction_candidates.append({
                "i": i,
                "score": score,
                "which": which,
                "mode": "junction",
                "lseq": lseq, "lpos": lpos, "llen": llen, "ltm": ltm, "l_ov": l_ov,
                "rseq": rseq, "rpos": rpos, "rlen": rlen, "rtm": rtm, "r_ov": r_ov,
                "ptm": ptm, "ppen": ppen, "psz": psz,
                "breakdown": breakdown,
                "kv": kv,
                "template": template,
                "jpos": jpos,
                "left_part_len": left_part_len,
                "map_mode": map_mode,
                "template_tx_start": template_start_tx,
                "tx_id": tx_id,
                "tx_version": tx_version,
                "eL_id": eL.exon_id,
                "eR_id": eR.exon_id,
                "eL_num": eL_num,
                "eR_num": eR_num,
                "total_junctions": total_junctions,
                "tx_len": tx_len,
                "junction_tx_pos": junction_tx_pos,
                "eL_len": eL_len,
                "eR_len": eR_len,
                "eL_tx_start": eL_tx_start,
                "eR_tx_start": eR_tx_start,
            })

        if len(junction_candidates) < int(args.min_candidates_per_junction):
            last_err = (
                f"此 junction（{eL.exon_id}->{eR.exon_id}）候選不足"
                f"（{len(junction_candidates)} < {args.min_candidates_per_junction}）。"
            )
            continue

        all_candidates.extend(junction_candidates)

    if not junction_mode:
        if args.span != "ANY" or int(args.min_junction_overlap) > 0:
            print("Warning: single-exon mode; --span/--min-junction-overlap will be ignored.", file=os.sys.stderr)

        if args.region:
            template = region_template
            template_start_tx = template_tx_start
        else:
            template, template_start_tx = build_center_template(cdna_seq, args.flank)

        kv = run_primer3_pairs(
            primer3_core=primer3_core,
            template=template,
            prod_min=args.prod_min,
            prod_max=args.prod_max,
            num_return=args.num_return,
        )

        single_candidates = []
        for i in range(args.num_return):
            parsed = parse_pair(kv, i)
            if not parsed:
                continue
            lseq, lpos, llen, ltm, rseq, rpos, rlen, rtm, ptm, ppen, psz = parsed

            score, breakdown = score_pair(
                kv=kv, i=i,
                lseq=lseq, lpos=lpos, llen=llen, ltm=ltm,
                rseq=rseq, rpos=rpos, rlen=rlen, rtm=rtm,
                pair_tm=ptm, product_size=psz,
                junction_pos=0,
                opt_tm=60.0,
                opt_size=args.opt_amplicon,
            )

            single_candidates.append({
                "i": i,
                "score": score,
                "which": "N/A",
                "mode": "single",
                "lseq": lseq, "lpos": lpos, "llen": llen, "ltm": ltm, "l_ov": 0,
                "rseq": rseq, "rpos": rpos, "rlen": rlen, "rtm": rtm, "r_ov": 0,
                "ptm": ptm, "ppen": ppen, "psz": psz,
                "breakdown": breakdown,
                "kv": kv,
                "template": template,
                "jpos": -1,
                "left_part_len": None,
                "map_mode": "offset",
                "template_tx_start": template_start_tx,
                "tx_id": tx_id,
                "tx_version": tx_version,
                "eL_id": None,
                "eR_id": None,
                "eL_num": None,
                "eR_num": None,
                "total_junctions": 0,
                "tx_len": tx_len,
                "junction_tx_pos": None,
                "eL_len": None,
                "eR_len": None,
                "eL_tx_start": None,
                "eR_tx_start": None,
            })

        if len(single_candidates) < int(args.min_candidates_per_junction):
            last_err = (
                f"single-exon mode 候選不足"
                f"（{len(single_candidates)} < {args.min_candidates_per_junction}）。"
            )
        else:
            all_candidates.extend(single_candidates)

    if not all_candidates:
        raise RuntimeError(last_err or "多個 junction 嘗試後仍無法找到符合條件的 primer pair。")

    all_candidates.sort(key=lambda x: x["score"], reverse=True)
    topn = all_candidates[: max(1, int(args.top))]

    lines: List[str] = []
    lines.append("Top candidates (ranked across junctions):")
    # Print shared transcript/junction info once (from the top-ranked candidate)
    c0 = topn[0]
    if user_sequence:
        lines.append("Input: user sequence")
    else:
        lines.append(f"Species: {species}")
        if gene_name:
            lines.append(f"Gene name: {gene_name}")
        lines.append(f"Gene: {gene}")
    if c0["tx_version"] is not None:
        lines.append(f"Transcript（挑選）: {c0['tx_id']}.{c0['tx_version']}")
    else:
        lines.append(f"Transcript（挑選）: {c0['tx_id']}")
    if c0["mode"] == "junction":
        lines.append(
            f"Junction（exon）: {c0['eL_num']} ({c0['eL_id']}) -> {c0['eR_num']} ({c0['eR_id']})"
        )
        lines.append(f"Total junctions in transcript: {c0['total_junctions']}")
        lines.append(f"Transcript 長度: {c0['tx_len']} bp（junction cDNA 位置: {c0['junction_tx_pos']}）")
    else:
        lines.append("Mode: single-exon (no junction)")
        lines.append("Exons in transcript: 1")
        lines.append(f"Transcript 長度: {c0['tx_len']} bp")
    if args.region:
        lines.append(f"Region (cDNA, 1-based): {args.region}")
    lines.append("")
    for rank, c in enumerate(topn, start=1):
        i = c["i"]
        lseq = c["lseq"]; lpos = c["lpos"]; llen = c["llen"]; ltm = c["ltm"]
        rseq = c["rseq"]; rpos = c["rpos"]; rlen = c["rlen"]; rtm = c["rtm"]
        ptm  = c["ptm"];  psz  = c["psz"];  which = c["which"]
        ppen = c["ppen"]
        kv = c["kv"]

        best = JunctionDesign(
            gene=gene,
            transcript_id=c["tx_id"],
            transcript_version=c["tx_version"],
            exon_left=c["eL_id"],
            exon_right=c["eR_id"],
            exon_left_num=c["eL_num"],
            exon_right_num=c["eR_num"],
            left_primer=lseq,
            right_primer=rseq,
            left_tm=ltm,
            right_tm=rtm,
            pair_tm=ptm,
            pair_penalty=ppen,
            product_size=psz,
            template=c["template"],
            junction_pos=c["jpos"],
            left_pos=lpos,
            left_len=llen,
            right_pos=rpos,
            right_len=rlen,
            which_spans=which,
            junction_overlap_left=c["l_ov"],
            junction_overlap_right=c["r_ov"],
            score=c["score"],
        )

        lines.append(f"== Rank {rank} | score={best.score:.3f} ==")
        if c["map_mode"] == "junction":
            eL_len = c["eL_len"]
            eR_len = c["eR_len"]
            eL_tx_start = c["eL_tx_start"]
            eR_tx_start = c["eR_tx_start"]
            left_part_len = c["left_part_len"]
            l_tx_start = map_template_pos_to_tx(
                best.left_pos,
                left_exon_len=eL_len,
                right_exon_len=eR_len,
                left_exon_tx_start=eL_tx_start,
                right_exon_tx_start=eR_tx_start,
                left_part_len=left_part_len,
            )
            l_tx_end = map_template_pos_to_tx(
                best.left_pos + best.left_len - 1,
                left_exon_len=eL_len,
                right_exon_len=eR_len,
                left_exon_tx_start=eL_tx_start,
                right_exon_tx_start=eR_tx_start,
                left_part_len=left_part_len,
            )
            r_tx_start = map_template_pos_to_tx(
                best.right_pos,
                left_exon_len=eL_len,
                right_exon_len=eR_len,
                left_exon_tx_start=eL_tx_start,
                right_exon_tx_start=eR_tx_start,
                left_part_len=left_part_len,
            )
            r_tx_end = map_template_pos_to_tx(
                best.right_pos + best.right_len - 1,
                left_exon_len=eL_len,
                right_exon_len=eR_len,
                left_exon_tx_start=eL_tx_start,
                right_exon_tx_start=eR_tx_start,
                left_part_len=left_part_len,
            )
        else:
            template_start_tx = c["template_tx_start"] or 0
            l_tx_start = map_template_pos_to_tx_offset(best.left_pos, template_start_tx)
            l_tx_end = map_template_pos_to_tx_offset(best.left_pos + best.left_len - 1, template_start_tx)
            r_tx_start = map_template_pos_to_tx_offset(best.right_pos, template_start_tx)
            r_tx_end = map_template_pos_to_tx_offset(best.right_pos + best.right_len - 1, template_start_tx)
        lines.append(
            f"LEFT : {best.left_primer}   (cDNA={l_tx_start+1}-{l_tx_end+1}, len={best.left_len}, Tm={best.left_tm:.2f})"
        )
        lines.append(
            f"RIGHT: {best.right_primer}  (cDNA={r_tx_start+1}-{r_tx_end+1}, len={best.right_len}, Tm={best.right_tm:.2f})"
        )
        l_gc = gc_percent(best.left_primer)
        r_gc = gc_percent(best.right_primer)
        lines.append(f"GC%  : LEFT {l_gc:.1f} | RIGHT {r_gc:.1f}")
        if c["mode"] == "junction":
            lines.append(
                f"Junction overlap (bp): LEFT {best.junction_overlap_left} | RIGHT {best.junction_overlap_right}  (min={args.min_junction_overlap}, rule={args.span})"
            )
        else:
            lines.append("Junction overlap (bp): N/A (single-exon mode)")

        l_self_any = pick_metric(
            kv,
            [f"PRIMER_LEFT_{i}_SELF_ANY_TH", f"PRIMER_LEFT_{i}_SELF_ANY"],
        )
        l_self_end = pick_metric(
            kv,
            [f"PRIMER_LEFT_{i}_SELF_END_TH", f"PRIMER_LEFT_{i}_SELF_END"],
        )
        l_hairpin = pick_metric(
            kv,
            [f"PRIMER_LEFT_{i}_HAIRPIN_TH", f"PRIMER_LEFT_{i}_HAIRPIN"],
        )
        r_self_any = pick_metric(
            kv,
            [f"PRIMER_RIGHT_{i}_SELF_ANY_TH", f"PRIMER_RIGHT_{i}_SELF_ANY"],
        )
        r_self_end = pick_metric(
            kv,
            [f"PRIMER_RIGHT_{i}_SELF_END_TH", f"PRIMER_RIGHT_{i}_SELF_END"],
        )
        r_hairpin = pick_metric(
            kv,
            [f"PRIMER_RIGHT_{i}_HAIRPIN_TH", f"PRIMER_RIGHT_{i}_HAIRPIN"],
        )
        p_compl_any = pick_metric(
            kv,
            [f"PRIMER_PAIR_{i}_COMPL_ANY_TH", f"PRIMER_PAIR_{i}_COMPL_ANY"],
        )
        p_compl_end = pick_metric(
            kv,
            [f"PRIMER_PAIR_{i}_COMPL_END_TH", f"PRIMER_PAIR_{i}_COMPL_END"],
        )

        if any([l_self_any, l_self_end, l_hairpin, r_self_any, r_self_end, r_hairpin, p_compl_any, p_compl_end]):
            if l_self_any or l_self_end or l_hairpin:
                lines.append(
                    "LEFT  dimer/hairpin: "
                    + ", ".join(
                        x
                        for x in [
                            f"self_any={l_self_any}" if l_self_any else "",
                            f"self_end={l_self_end}" if l_self_end else "",
                            f"hairpin={l_hairpin}" if l_hairpin else "",
                        ]
                        if x
                    )
                )
            if r_self_any or r_self_end or r_hairpin:
                lines.append(
                    "RIGHT dimer/hairpin: "
                    + ", ".join(
                        x
                        for x in [
                            f"self_any={r_self_any}" if r_self_any else "",
                            f"self_end={r_self_end}" if r_self_end else "",
                            f"hairpin={r_hairpin}" if r_hairpin else "",
                        ]
                        if x
                    )
                )
            if p_compl_any or p_compl_end:
                lines.append(
                    "PAIR  heterodimer: "
                    + ", ".join(
                        x
                        for x in [
                            f"compl_any={p_compl_any}" if p_compl_any else "",
                            f"compl_end={p_compl_end}" if p_compl_end else "",
                        ]
                        if x
                    )
                )

        if best.pair_tm is not None:
            lines.append(f"Pair product Tm: {best.pair_tm:.2f}")
        if best.pair_penalty is not None:
            lines.append(f"Pair penalty (primer3): {best.pair_penalty:.3f}")
        lines.append(f"預期產物長度: {best.product_size} bp")
        amp_start = best.left_pos
        amp_end = best.right_pos + 1
        if 0 <= amp_start < amp_end <= len(best.template):
            amplicon = best.template[amp_start:amp_end]
            lines.append(f"預期產物序列: {amplicon}")
        if c["mode"] == "junction":
            lines.append(f"跨 junction 的引子: {best.which_spans}")
        else:
            lines.append("跨 junction 的引子: N/A")
        lines.append("")

        if c["mode"] == "junction" and 0 <= best.junction_pos <= len(best.template):
            lo = max(0, best.junction_pos - 25)
            hi = min(len(best.template), best.junction_pos + 25)
            ctx = best.template[lo:hi]
            marker = " " * (best.junction_pos - lo) + "^ junction"
            lines.append("Junction context（左右各 25bp）:")
            lines.append(ctx)
            lines.append(marker)

    output = "\n".join(lines)
    print(output)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(output + "\n")
    return

if __name__ == "__main__":
    main()
