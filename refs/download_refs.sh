#!/bin/bash
set -e

REF_DIR="/mnt/scratch/DM/cfdna/refs"

cd $REF_DIR

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes

echo "Downloaded to $REF_DIR!"
