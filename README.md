qp_primer_design.py

Purpose
- Design RT-qPCR primers (junction or single-exon) using primer3_core.
- Optional BLAST-to-Ensembl resolution for sequence-only input.
- Optional QC using seqkit amplicon against transcriptome FASTA.

Requirements
- Python 3
- primer3_core (compiled)
- seqkit (only if --qc is used)

Installation (quick)
- Ensure primer3_core is compiled and executable.
- Install seqkit only if you plan to use QC.

Basic usage (gene mode)
python qp_primer_design.py \
  --gene TP53 \
  --species homo_sapiens \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt

Sequence mode (U will be converted to T)
python qp_primer_design.py \
  --sequence ACTG... \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt

BLAST-assisted sequence mode (maps to Ensembl when possible)
python qp_primer_design.py \
  --sequence ACTG... \
  --blast \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt

Region-restricted design (cDNA 1-based inclusive)
python qp_primer_design.py \
  --gene TP53 \
  --species homo_sapiens \
  --region 200-400 \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt

Single-exon mode
- If a transcript has only one exon, the script automatically switches to single-exon mode.
- You can also force allow with --allow-single-exon.

QC (canonical unique policy)
python qp_primer_design.py \
  --gene TP53 \
  --species homo_sapiens \
  --qc \
  --primer3 /abs/path/to/primer3_core \
  --out /abs/path/output.txt

QC reference behavior
- If --qc-ref-fasta is provided, that file is used.
- Otherwise, the script auto-downloads Ensembl cdna.all.fa.gz for the given species.
- The reference is downloaded into the same folder as --out and cached for reuse.

Common options
--gene <symbol>
--species <ensembl_species>
--sequence <DNA>
--blast
--primer3 <path>
--region <start-end>
--allow-single-exon
--qc
--qc-ref-fasta <path>
--ensembl-release <int or current>
--out <path>
--top <int>

Outputs
- Main report: --out
- QC report (if --qc): <out>.qc.tsv

Design modes
- Junction mode (default): If transcript has >= 2 exons, primers are designed on exon-exon junctions.
- Single-exon mode: If transcript has only 1 exon, it switches automatically; junction constraints are ignored.
- Sequence-only mode: Provide --sequence; transcript/gene fields are not used.

QC policy (canonical_unique)
- PASS only if:
  - Exactly 1 amplicon is found in the transcriptome reference.
  - The amplicon is on the canonical transcript.
- Otherwise FAIL and the next candidate is tried.

3' end filtering rules (hard filter)
Stage 1 (last 3 bp)
- 3' terminal base must be G/C.
- GC count in last 3 bp >= 1.
- Last 3 bp cannot be GGG or CCC.
- No 3' homopolymer run >= 4.

Stage 2 (last 5 bp)
- GC count in last 5 bp between 1 and 3.
- Max run (any nt) < 4.
- Max run (G/C) < 3.
- No strong 3' self-complementarity (simple check).
- No strong 3' cross-dimer between primers (simple check).

BLAST behavior (sequence mode with --blast)
- Uses NCBI BLAST (db=nt).
- Tries top N hits (default 5) to map to Ensembl gene ID.
- If no Ensembl mapping, falls back to sequence-only mode and prints a warning.

Notes
- Region coordinates use transcript cDNA, 1-based inclusive (e.g., --region 200-400).
- If QC is enabled and reference is not provided, the script auto-downloads Ensembl cdna.all.fa.gz into the output directory and reuses it if already present.
