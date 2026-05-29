#!/bin/bash
# RUN ON HOST (sambamba/samtools live on the host, not in the analysis container).
# Per-base cfDNA coverage in each tissue's super-enhancer +-1 kb windows, collapsed to a
# site-relative OFFSET PROFILE (-1000..+1000 from the SE center) -- the aggregated coverage
# "wiggle". Library size (total mapped reads) is recorded for read-scaling; the read-scaling
# and the cancer / matched-normal per-offset division happen in 02_plot_coverage.R.
#
#   for each sample x tissue in coverage_samplesheet.csv:
#     sambamba depth base (MAPQ>=30, no dup/secondary/qcfail/unmapped)
#       | aggregate_offsets.py SITES SAMPLE TISSUE TOTAL_READS
#       -> coverage/{sample}.{tissue}.offset_cov.tsv   (offset, sum_cov, n_covered)
#
# The 'tissues' column restricts which site sets a sample is measured against ('all' = every
# tissue; ';'-separated subset otherwise). To TEST one tumor-normal pair first, set the env
# vars ONLY_SAMPLES and ONLY_TISSUES (comma/space lists), e.g.:
#     ONLY_SAMPLES="IC15 IH02" ONLY_TISSUES=lung bash 01_coverage.sh
# Idempotent: existing outputs are skipped.
set -euo pipefail

ROOT="${ROOT:-/mnt/scratch/DM/cfdna}"
COO="$ROOT/analysis/cell_of_origin"
# Preprocessed BAMs from the finale pipeline: coordinate-sorted, MarkDuplicates-flagged, indexed.
WORKDIR="$ROOT/cfdna-finale-snakemake/work"
SHEET="$COO/coverage_samplesheet.csv"
SITESDIR="$COO/sites"
OUTDIR="$COO/coverage"; mkdir -p "$OUTDIR"
AGG="$COO/aggregate_offsets.py"
SAMBAMBA="${SAMBAMBA:-sambamba}"; SAMTOOLS="${SAMTOOLS:-samtools}"
THREADS="${THREADS:-8}"
FILTER="mapping_quality >= 30 and not duplicate and not secondary_alignment and not failed_quality_control and not unmapped"

# Optional test filters (comma/space separated); empty = use the sheet as-is.
ONLY_SAMPLES="${ONLY_SAMPLES:-}"; ONLY_TISSUES="${ONLY_TISSUES:-}"

command -v "$SAMBAMBA" >/dev/null || { echo "ERROR: sambamba not found (set \$SAMBAMBA)"; exit 1; }
command -v "$SAMTOOLS" >/dev/null || { echo "ERROR: samtools not found (set \$SAMTOOLS)"; exit 1; }
[[ -s "$SHEET" ]] || { echo "ERROR: control sheet not found: $SHEET"; exit 1; }
[[ -s "$AGG" ]]   || { echo "ERROR: aggregator not found: $AGG"; exit 1; }

mapfile -t ALL_TISSUES < <(ls "$SITESDIR"/*.bed 2>/dev/null | xargs -n1 basename | sed 's/\.bed$//')
[[ ${#ALL_TISSUES[@]} -gt 0 ]] || { echo "ERROR: no sites/*.bed -- run 00_build_sites.sh first"; exit 1; }

want=" ${ONLY_SAMPLES//,/ } "   # padded so " IC15 " word-matches

# control sheet columns: sample_id,sample_group,library_type,coverage,tissues
tail -n +2 "$SHEET" | tr -d '\r' | while IFS=, read -r sample group lib cov tissues _rest; do
  [[ -n "$sample" ]] || continue
  [[ -z "$ONLY_SAMPLES" || "$want" == *" $sample "* ]] || continue
  bam="$WORKDIR/$sample/$sample.md.bam"
  [[ -s "$bam" ]] || { echo "  ! missing BAM $bam -- skip $sample"; continue; }

  reads=$("$SAMTOOLS" idxstats "$bam" | awk '{m+=$3} END{print m+0}')   # library size

  if [[ -n "$ONLY_TISSUES" ]]; then
    IFS=', ' read -ra sel <<< "$ONLY_TISSUES"
  elif [[ "${tissues,,}" == "all" || -z "$tissues" ]]; then
    sel=("${ALL_TISSUES[@]}")
  else
    IFS=';' read -ra sel <<< "$tissues"
  fi

  for tissue in "${sel[@]}"; do
    [[ -n "$tissue" ]] || continue
    bed="$SITESDIR/$tissue.bed"
    [[ -s "$bed" ]] || { echo "  ! no $bed -- skip $sample/$tissue"; continue; }
    out="$OUTDIR/$sample.$tissue.offset_cov.tsv"
    [[ -s "$out" ]] && { echo "  = $sample/$tissue exists"; continue; }
    echo "  + $sample/$tissue (reads=$reads)"
    "$SAMBAMBA" depth base -t "$THREADS" -F "$FILTER" -L "$bed" "$bam" \
      | python3 "$AGG" "$bed" "$sample" "$tissue" "$reads" > "$out.tmp" \
      && mv "$out.tmp" "$out"
  done
done
echo "offset profiles -> $OUTDIR"
