#!/bin/bash
# Build per-tissue TFBS site sets for the cell-of-origin coverage analysis.
#
#   GTRD TF-binding-site centers (hg38)  INTERSECT  tissue-specific super-enhancers (SEdb, hg38)
#     -> isolate (drop sites whose ±1 kb windows would overlap)  -> liftOver hg38->hg19
#     -> strip 'chr', restrict to 1-22/X/Y, expand to ±1 kb windows, bound-check
#   = sites/{tissue}.bed  (chrom, start, end, center)  used by 01_coverage.sh
#
# Tissue -> SEdb tissue-type mapping is in tissue_map.tsv. SEs are taken from normal-lineage
# samples (Biosample type = Tissue / Primary cell), not cancer cell lines. Requires the
# references fetched by refs/download_refs_coo.sh. Runs in-container (bedtools + liftOver).
set -euo pipefail

ROOT="${ROOT:-/mnt/scratch/DM/cfdna}"
COO="$ROOT/analysis/cell_of_origin"
REF="$ROOT/refs/coo"
GTRD="$REF/gtrd_hg38_tf_centers.bed.gz"
SEINFO="$REF/sedb3_human_sample_information.txt"
LIFTOVER="$REF/liftOver"
CHAIN="$REF/hg38ToHg19.over.chain.gz"
SIZES="$ROOT/refs/hg19.chrom.sizes"
SE_URL_BASE="http://www.licpathway.net/sedb/download_v3/SE_hg38_bed"
TMAP="$COO/tissue_map.tsv"

SITES="$COO/sites"; SETMP="$COO/se_tmp"; SEDL="$SETMP/se_dl"
mkdir -p "$SITES" "$SETMP" "$SEDL"

WINDOW=1000; MINGAP=2000; MAXSAMPLES="${MAXSAMPLES:-60}"; MAXSITES="${MAXSITES:-20000}"

for f in "$GTRD" "$SEINFO" "$LIFTOVER" "$CHAIN" "$SIZES" "$TMAP"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

keep_re='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'
manifest="$SITES/manifest.tsv"
printf "tissue\tn_samples\tn_se_union\tn_centers_hg38\tn_sites_final\n" > "$manifest"

tail -n +2 "$TMAP" | while IFS=$'\t' read -r tissue regex note; do
  [[ -n "$tissue" ]] || continue
  echo "=== $tissue ($regex) ==="

  # 1) select normal-lineage SEdb samples for this tissue
  mapfile -t ids < <(awk -F'\t' -v re="$regex" \
      'NR>1 && tolower($5) ~ re && tolower($4) ~ /tissue|primary/ {print $1}' "$SEINFO" \
      | sort -u | head -n "$MAXSAMPLES")
  echo "  samples: ${#ids[@]}"

  # 2) download each sample's SE bed (cached), pool chrom/start/end (skip header line)
  raw="$SETMP/$tissue.hg38.se.raw.bed"; : > "$raw"
  for id in "${ids[@]}"; do
    bed="$SEDL/${id}_SE_hg38.bed"
    if [[ ! -s "$bed" ]]; then
      curl -sL --max-time 60 -o "$bed" "$SE_URL_BASE/${id}_SE_hg38.bed" || true
    fi
    # valid SEdb bed begins with the se_chr header
    head -1 "$bed" 2>/dev/null | grep -q '^se_chr' \
      && tail -n +2 "$bed" | cut -f1-3 >> "$raw" || echo "    ! bad/empty $id"
  done

  # 3) SE union (hg38), standard chroms only
  seu="$SETMP/$tissue.hg38.se.bed"
  awk -v re="$keep_re" '$1 ~ re' "$raw" | sort -k1,1 -k2,2n | bedtools merge > "$seu"
  n_seu=$(wc -l < "$seu")

  # 4) GTRD centers within tissue SEs (hg38) -> point bed with index name
  cen38="$SETMP/$tissue.hg38.centers.bed"
  bedtools intersect -a <(zcat "$GTRD") -b "$seu" -u \
    | awk -v re="$keep_re" '$1 ~ re {print $1"\t"$2"\t"$3"\t"NR}' > "$cen38"
  n_cen=$(wc -l < "$cen38")

  # 5) liftOver centers hg38 -> hg19
  cen19="$SETMP/$tissue.hg19.centers.bed"; unmap="$SETMP/$tissue.unmapped.bed"
  "$LIFTOVER" "$cen38" "$CHAIN" "$cen19" "$unmap" >/dev/null 2>&1 || true

  # 6) hg19: strip chr, isolate (greedy >= MINGAP), cap, expand ±WINDOW, bound-check
  iso="$SETMP/$tissue.hg19.iso.txt"
  awk -v re="$keep_re" '$1 ~ re {sub(/^chr/,"",$1); print $1"\t"$2}' "$cen19" \
    | sort -k1,1 -k2,2n \
    | awk -v g="$MINGAP" 'BEGIN{pc="";pp=-1e9} {if($1!=pc){pc=$1;pp=-1e9} if($2-pp>=g){print $1"\t"$2; pp=$2}}' \
    > "$iso"
  n_iso=$(wc -l < "$iso")

  step=1
  if (( n_iso > MAXSITES )); then step=$(( (n_iso + MAXSITES - 1) / MAXSITES )); fi

  awk -v w="$WINDOW" -v step="$step" -v sizesf="$SIZES" '
    BEGIN{ while((getline l < sizesf)>0){ split(l,a,"\t"); sub(/^chr/,"",a[1]); sz[a[1]]=a[2] } }
    { if((NR-1)%step!=0) next
      c=$2; s=c-w; e=c+w+1
      if(s>=0 && (($1 in sz)? e<=sz[$1] : 1)) print $1"\t"s"\t"e"\t"c }
  ' "$iso" | sort -k1,1 -k2,2n > "$SITES/$tissue.bed"
  n_final=$(wc -l < "$SITES/$tissue.bed")

  echo "  SE union=$n_seu  centers(hg38)=$n_cen  isolated=$n_iso  final sites=$n_final (step=$step)"
  printf "%s\t%d\t%d\t%d\t%d\n" "$tissue" "${#ids[@]}" "$n_seu" "$n_cen" "$n_final" >> "$manifest"
done

echo; echo "=== manifest ==="; column -t "$manifest"
echo "sites -> $SITES"
