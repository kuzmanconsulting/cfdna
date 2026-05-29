#!/bin/bash
# Reference data for the cell-of-origin cfDNA coverage analysis. RUN ON HOST (good bandwidth).
#
#   1. SEdb 3.0 human sample-information table (sample -> tissue / biosample type).
#   2. The tissues our cfDNA panel needs, derived by joining:
#        coverage_samplesheet.csv  ->  snyder2016 disease  ->  tumor_types.csv tissue
#   3. Per tissue, the NORMAL-lineage SEdb samples (Biosample type = Tissue or Primary cell;
#      cancer cell lines excluded) -> se_samples.tsv  (sample_id, tissue).
#   4. Each selected sample's whole super-enhancer BED (hg19), cached in se_hg19/.
#
# We fetch SE (super-enhancer) BEDs, not SE_ele: coverage is profiled at the super-enhancer
# center +-1 kb, so only the SE coordinates are needed. BAMs are GRCh37/hg19 and SEdb serves
# hg19, so there is no liftOver step. 00_build_sites.sh consumes se_samples.tsv + se_hg19/.
set -euo pipefail

ROOT="${ROOT:-/mnt/scratch/DM/cfdna}"
REF_DIR="$ROOT/cell-of-origin/refs"
SEDB_BASE="${SEDB_BASE:-http://www.licpathway.net/sedb/download_v3}"
SEINFO="$REF_DIR/sedb3_human_sample_information.txt"
SHEET="$ROOT/cell-of-origin/coverage_samplesheet.csv"
META="$ROOT/snyder2016_metadata_GSE.csv"
TUMOR="$ROOT/analysis/tumor_types.csv"
SAMPLES="$REF_DIR/se_samples.tsv"
SEDIR="$REF_DIR/se_hg19"
mkdir -p "$SEDIR"

# tumor_types tissue -> SEdb 'Tissue type' regex. Default is the tissue name as a substring;
# only breast/liver need a wider pattern to catch SEdb's vocabulary (mammary, hepatocyte).
declare -A TIS_RE=( [breast]='breast|mammary' [liver]='liver|hepat' )

for f in "$SHEET" "$META" "$TUMOR"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

# 1. SEdb sample table (downloaded if absent)
[[ -s "$SEINFO" ]] || \
  curl -sL --connect-timeout 20 --max-time 120 -o "$SEINFO" "$SEDB_BASE/Human_sample_information_sedb3.txt"
echo "[refs] SEdb samples: $(($(wc -l < "$SEINFO") - 1))"

# 2. panel tissues: samplesheet sample -> metadata disease -> tumor_types tissue
mapfile -t TISSUES < <(awk -F',' '
  FILENAME==ARGV[1] && FNR>1 { sub(/\r$/,"",$3); dt[$1]=$3 }                 # tumor_types: disease->tissue
  FILENAME==ARGV[2] && FNR>1 { sub(/\r$/,"",$6); sd[$1]=$6 }                 # metadata:    sample ->disease
  FILENAME==ARGV[3] && FNR>1 { sub(/\r$/,"",$1); if(dt[sd[$1]]!="") print dt[sd[$1]] }   # samplesheet sample
' "$TUMOR" "$META" "$SHEET" | sort -u)
echo "[refs] panel tissues: ${TISSUES[*]}"

# 3. normal-lineage SEdb samples per tissue (cols: ID=1, Biosample type=4, Tissue type=5)
: > "$SAMPLES"
for t in "${TISSUES[@]}"; do
  awk -F'\t' -v re="${TIS_RE[$t]:-$t}" -v t="$t" '
    NR>1 && tolower($4) ~ /^tissue$|^primary cell$/ && tolower($5) ~ re { print $1"\t"t }
  ' "$SEINFO" >> "$SAMPLES"
done
sort -u -o "$SAMPLES" "$SAMPLES"
echo "[refs] reference samples per tissue:"; cut -f2 "$SAMPLES" | sort | uniq -c

# 4. fetch each sample's super-enhancer BED (cached; re-runs skip existing files)
while IFS=$'\t' read -r id _tissue; do
  bed="$SEDIR/${id}_SE_hg19.bed"
  [[ -s "$bed" ]] && continue
  curl -sL --connect-timeout 20 --max-time 120 -o "$bed" "$SEDB_BASE/SE_hg19_bed/${id}_SE_hg19.bed" || true
  head -1 "$bed" 2>/dev/null | grep -q '^se_chr' || { echo "  ! bad/empty $id"; rm -f "$bed"; }
done < "$SAMPLES"

echo "[refs] done -> $SEDIR ($(ls "$SEDIR" | wc -l) SE BEDs cached); samples -> $SAMPLES"
