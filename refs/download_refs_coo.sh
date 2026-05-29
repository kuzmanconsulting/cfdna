#!/bin/bash
# Reference data for the cell-of-origin cfDNA coverage analysis.
# Run on host (good bandwidth); writes to the shared refs/coo mount.
#
#   GTRD 19.10 human TF meta-clusters (hg38) -> slimmed to TF-binding-site centers
#   SEdb 3.0 human sample-information table (sample -> tissue map)
#   UCSC liftOver binary + hg38->hg19 chain (final sites are lifted to hg19 to match BAMs)
#
# BAMs are GRCh37/hg19; GTRD is hg38, so intersection happens in hg38 and the final
# per-tissue site sets are lifted to hg19 by 00_build_sites.sh.
set -euo pipefail

REF_DIR="/mnt/scratch/DM/cfdna/refs/coo"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

GTRD_URL="https://gtrd.biouml.org/downloads/19.10/chip-seq/Homo%20sapiens_meta_clusters.interval.gz"
SEDB_INFO_URL="http://www.licpathway.net/sedb/download_v3/Human_sample_information_sedb3.txt"
LIFTOVER_URL="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"
CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

# meta_clusters columns: CHROM START END summit uniprotId tfTitle cell.set treatment.set
#                        exp.set peak-caller.set peak-caller.count exp.count peak.count
# Keep meta-clusters supported by >=2 experiments; emit 1 bp summit center (START+summit)
# with the TF title. ~4.5 GB download streamed, never stored uncompressed.
# Re-fetch unless a *complete* gzip already exists (gzip -t catches truncated downloads).
# Write to .tmp and rename, so an interrupted run never leaves a bad final file.
if ! gzip -t gtrd_hg38_tf_centers.bed.gz 2>/dev/null; then
  echo "[refs] streaming + slimming GTRD meta-clusters (hg38) ..."
  curl -sL "$GTRD_URL" \
    | zcat \
    | awk -F'\t' 'NR>1 && $12>=2 { s=$4; c=(s>=0)?($2+s):int(($2+$3)/2); print $1"\t"c"\t"c+1"\t"$6 }' \
    | gzip > gtrd_hg38_tf_centers.bed.gz.tmp \
    && mv gtrd_hg38_tf_centers.bed.gz.tmp gtrd_hg38_tf_centers.bed.gz
fi
echo "[refs] GTRD centers: $(zcat gtrd_hg38_tf_centers.bed.gz | wc -l) sites"

[[ -s sedb3_human_sample_information.txt ]] || curl -sL -o sedb3_human_sample_information.txt "$SEDB_INFO_URL"
echo "[refs] SEdb samples: $(($(wc -l < sedb3_human_sample_information.txt) - 1))"

[[ -x liftOver ]] || { curl -sL -o liftOver "$LIFTOVER_URL"; chmod +x liftOver; }
[[ -s hg38ToHg19.over.chain.gz ]] || curl -sL -o hg38ToHg19.over.chain.gz "$CHAIN_URL"

echo "[refs] done -> $REF_DIR"
