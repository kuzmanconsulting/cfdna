#!/bin/bash
# Build per-tissue super-enhancer site sets for the cell-of-origin coverage analysis.
#
#   SEdb 3.0 super-enhancers (hg19) for each tissue's reference samples (se_samples.tsv)
#     -> center = midpoint of each super-enhancer
#     -> pool over the tissue's samples, keep chr1-22/X/Y, strip 'chr' (BAMs are b37)
#     -> isolate (drop centers whose +-1 kb windows would overlap)
#     -> expand to +-1 kb windows, bound-check against chrom sizes
#   = sites/{tissue}.bed  (chrom, start, end, center)  consumed by 01_coverage.sh
#
# Reference SE BEDs and se_samples.tsv come from refs/download_refs_coo.sh (run that first).
# hg19-native: no liftOver, no bedtools -- awk only.
set -euo pipefail

ROOT="${ROOT:-/mnt/scratch/DM/cfdna}"
COO="$ROOT/analysis/cell_of_origin"
REF="$ROOT/refs/coo"
SAMPLES="$REF/se_samples.tsv"
SIZES="$ROOT/refs/hg19.chrom.sizes"
SEDIR="$REF/se_hg19"
SITES="$COO/sites"; mkdir -p "$SITES"

WINDOW=1000; MINGAP=2000

for f in "$SAMPLES" "$SIZES"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f (run download_refs_coo.sh first)"; exit 1; }
done

keep_re='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'
manifest="$SITES/manifest.tsv"
printf "tissue\tn_samples\tn_superenhancers\tn_sites\n" > "$manifest"

for tissue in $(cut -f2 "$SAMPLES" | sort -u); do
  echo "=== $tissue ==="
  raw=$(mktemp); ns=0
  while IFS=$'\t' read -r id t; do
    [[ "$t" == "$tissue" ]] || continue
    bed="$SEDIR/${id}_SE_hg19.bed"
    [[ -s "$bed" ]] || { echo "  ! missing $bed"; continue; }
    awk -F'\t' 'NR>1 {print $1"\t"int(($2+$3)/2)}' "$bed" >> "$raw"   # se_chr, midpoint
    ns=$((ns+1))
  done < "$SAMPLES"
  n_se=$(wc -l < "$raw")

  # standard chroms -> strip chr -> sort/unique -> isolate (>= MINGAP so +-1 kb windows
  # never overlap) -> expand +-WINDOW -> bound-check against chrom sizes (bare contigs)
  awk -v re="$keep_re" '$1 ~ re {sub(/^chr/,"",$1); print $1"\t"$2}' "$raw" \
    | sort -k1,1 -k2,2n -u \
    | awk -v g="$MINGAP" 'BEGIN{pc="";pp=-1e9}{if($1!=pc){pc=$1;pp=-1e9} if($2-pp>=g){print;pp=$2}}' \
    | awk -v w="$WINDOW" -v sizesf="$SIZES" '
        BEGIN{ while((getline l < sizesf)>0){split(l,a,"\t"); sub(/^chr/,"",a[1]); sz[a[1]]=a[2]} }
        { s=$2-w; e=$2+w+1; if(s>=0 && (($1 in sz)? e<=sz[$1] : 1)) print $1"\t"s"\t"e"\t"$2 }' \
    | sort -k1,1 -k2,2n > "$SITES/$tissue.bed"
  n_sites=$(wc -l < "$SITES/$tissue.bed")
  rm -f "$raw"

  echo "  samples=$ns  super-enhancers=$n_se  sites=$n_sites"
  printf "%s\t%d\t%d\t%d\n" "$tissue" "$ns" "$n_se" "$n_sites" >> "$manifest"
done

echo; echo "=== manifest ==="; column -t "$manifest"
echo "sites -> $SITES"
