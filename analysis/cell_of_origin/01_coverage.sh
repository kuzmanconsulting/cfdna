#!/bin/bash
# RUN ON HOST. Per-base cfDNA coverage in the ±1 kb windows of each tissue's TFBS set,
# collapsed to a site-relative offset profile. sambamba is NOT in the analysis container,
# so this step runs on the host where sambamba is installed.
#
#   for each sample (samplesheet.csv) x tissue (sites/*.bed):
#     sambamba depth base (MAPQ>=30, no dup/secondary/qcfail/unmapped, all fragment sizes)
#       | aggregate_offsets.py  ->  coverage/{sample}.{tissue}.offset_cov.tsv
#
# Idempotent: existing outputs are skipped.
set -euo pipefail

ROOT="${ROOT:-/mnt/scratch/DM/cfdna}"
COO="$ROOT/analysis/cell_of_origin"
# Preprocessed BAMs from the finale pipeline: coordinate-sorted, MarkDuplicates-flagged,
# indexed. Required so the `not duplicate` filter and sambamba region access both work.
WORKDIR="$ROOT/cfdna-finale-snakemake/work"
SAMPLESHEET="$ROOT/cfdna-finale-snakemake/samplesheet.csv"
SITESDIR="$COO/sites"
OUTDIR="$COO/coverage"
AGG="$COO/aggregate_offsets.py"
SAMBAMBA="${SAMBAMBA:-sambamba}"
THREADS="${THREADS:-8}"
FILTER="mapping_quality >= 30 and not duplicate and not secondary_alignment and not failed_quality_control and not unmapped"

command -v "$SAMBAMBA" >/dev/null || { echo "ERROR: sambamba not found (set \$SAMBAMBA)"; exit 1; }
mkdir -p "$OUTDIR"

mapfile -t TISSUES < <(ls "$SITESDIR"/*.bed 2>/dev/null | xargs -n1 basename | sed 's/\.bed$//')
[[ ${#TISSUES[@]} -gt 0 ]] || { echo "ERROR: no sites/*.bed — run 00_build_sites.sh first"; exit 1; }
echo "tissues: ${TISSUES[*]}"

tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample srr _rest; do
  [[ -n "$sample" ]] || continue
  bam="$WORKDIR/$sample/$sample.md.bam"
  if [[ ! -s "$bam" ]]; then echo "  ! missing BAM $bam — skip $sample"; continue; fi
  for tissue in "${TISSUES[@]}"; do
    out="$OUTDIR/$sample.$tissue.offset_cov.tsv"
    [[ -s "$out" ]] && { echo "  = $sample/$tissue exists"; continue; }
    echo "  + $sample/$tissue"
    "$SAMBAMBA" depth base -t "$THREADS" -F "$FILTER" -L "$SITESDIR/$tissue.bed" "$bam" \
      | python3 "$AGG" "$SITESDIR/$tissue.bed" "$sample" "$tissue" > "$out.tmp" \
      && mv "$out.tmp" "$out"
  done
done
echo "coverage profiles -> $OUTDIR"
