# qp_primer_design.py

Unified RT-qPCR primer design pipeline with junction mode, single-exon mode, sequence mode, optional BLAST-to-Ensembl resolution, and optional QC via seqkit.

## 1. Install

### 1.1 Dependencies
- Python 3
- `primer3_core` (compiled and executable)
- `seqkit` (only if `--qc` is used)

### 1.2 Python packages

```sh
pip install requests
```

## 2. Quick Start

### 2.1 Gene mode

```sh
python qp_primer_design.py \
  --gene TP53 \
  --species homo_sapiens \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt
```

### 2.2 Sequence mode (U→T auto-converted)

```sh
python qp_primer_design.py \
  --sequence ACTG... \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt
```

### 2.3 BLAST-assisted sequence mode

```sh
python qp_primer_design.py \
  --sequence ACTG... \
  --blast \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt
```

### 2.4 Region-restricted design (cDNA 1-based inclusive)

```sh
python qp_primer_design.py \
  --gene TP53 \
  --species homo_sapiens \
  --region 200-400 \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt
```

## 3. Modes

- Junction mode (default): If transcript has >= 2 exons, primers are designed on exon–exon junctions.
- Single-exon mode: If transcript has only 1 exon, it switches automatically; junction constraints are ignored.
- Sequence-only mode: Provide `--sequence`; transcript/gene fields are not used.

## 4. QC (Canonical Unique)

Enable QC with `--qc`. This runs `seqkit amplicon` on a transcriptome FASTA and requires:
- Exactly 1 amplicon
- The amplicon is on the canonical transcript

### 4.1 QC reference
- If `--qc-ref-fasta` is provided, that file is used.
- Otherwise, the script auto-downloads Ensembl `cdna.all.fa.gz` for the given species and caches it in the output directory.

### 4.2 QC output
- Main report: `--out`
- QC report (if `--qc`): `<out>.qc.tsv`

## 5. 3' End Filtering (Hard Filter)

Stage 1 (last 3 bp):
- 3' terminal base must be G/C.
- GC count in last 3 bp >= 1.
- Last 3 bp cannot be GGG or CCC.
- No 3' homopolymer run >= 4.

Stage 2 (last 5 bp):
- GC count in last 5 bp between 1 and 3.
- Max run (any nt) < 4.
- Max run (G/C) < 3.
- No strong 3' self-complementarity (simple check).
- No strong 3' cross-dimer between primers (simple check).

## 6. Common Options

- `--gene <symbol>`: Ensembl gene symbol (used with `--species`).
- `--species <ensembl_species>`: Ensembl species name, e.g. `homo_sapiens`, `zea_mays`.
- `--sequence <seq>`: User-provided sequence (RNA or DNA). RNA will be auto-converted U→T.
- `--blast`: BLAST the sequence against NCBI `nt` and try to map to Ensembl.
- `--primer3 <path>`: Path to `primer3_core`.
- `--region <start-end>`: cDNA region (1-based, inclusive) to design within.
- `--allow-single-exon`: Allow single-exon mode (no junction).
- `--qc`: Enable QC with `seqkit amplicon`.
- `--qc-ref-fasta <path>`: Transcriptome FASTA for QC (cdna.all.fa.gz).
- `--ensembl-release <int or current>`: Ensembl release for auto-download (default `current`).
- `--out <path>`: Output report path.
- `--top <int>`: Number of top candidates to report.

## 7. Notes

- Region coordinates use transcript cDNA, 1-based inclusive (e.g., `--region 200-400`).
- BLAST uses NCBI `nt` and tries top N hits (default 5) to map to Ensembl.
- If BLAST→Ensembl mapping fails, it falls back to sequence-only mode.
