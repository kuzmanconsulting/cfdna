#!/usr/bin/env python3
"""Collapse `sambamba depth base` output into a site-relative coverage profile (the "wiggle").

Reads a sites BED (chrom, start, end, center) whose ±1 kb windows are isolated
(non-overlapping, super-enhancer centers >= 2 kb apart), so every genomic position maps to
exactly one center. Streams per-base depth on stdin and accumulates, per offset
(-1000..+1000 from the center), the summed coverage and the number of windows with a covered
base there. Aligning every super-enhancer on its center and averaging across sites yields the
aggregated coverage profile: mean coverage per offset = sum_cov / n_sites (the denominator is
the total window count, since zero-coverage bases are simply absent from sambamba's stream).

total_reads (library size) is recorded in the header so 02_plot_coverage.R can read-scale the
profile and divide each cancer by its matched normal, offset by offset.

Usage:  sambamba depth base ... <bam> | aggregate_offsets.py SITES.bed SAMPLE TISSUE TOTAL_READS
"""
import sys
from bisect import bisect_left

HALF = 1000

def main():
    sites_bed, sample, tissue, total_reads = sys.argv[1:5]

    centers = {}  # chrom -> sorted list of center positions
    n_sites = 0
    with open(sites_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, center = f[0], int(f[3])
            centers.setdefault(chrom, []).append(center)
            n_sites += 1
    for chrom in centers:
        centers[chrom].sort()

    sum_cov = [0] * (2 * HALF + 1)   # index = offset + HALF
    n_cov = [0] * (2 * HALF + 1)

    cur_chrom = None
    arr = []
    for line in sys.stdin:
        f = line.split("\t", 3)
        if len(f) < 3:
            continue
        chrom = f[0]
        try:
            pos = int(f[1]); cov = int(f[2])
        except ValueError:
            continue  # header / malformed
        if chrom != cur_chrom:
            arr = centers.get(chrom)
            cur_chrom = chrom
            if arr is None:
                continue
        if arr is None:
            continue
        i = bisect_left(arr, pos)
        # nearest center is arr[i-1] or arr[i]; windows are isolated so <=1 matches
        best = None
        if i < len(arr) and arr[i] - pos <= HALF:
            best = arr[i]
        if best is None and i > 0 and pos - arr[i - 1] <= HALF:
            best = arr[i - 1]
        if best is None:
            continue
        off = pos - best
        if -HALF <= off <= HALF:
            j = off + HALF
            sum_cov[j] += cov
            n_cov[j] += 1

    out = sys.stdout
    out.write(f"# n_sites={n_sites}\tsample={sample}\ttissue={tissue}\ttotal_reads={total_reads}\n")
    out.write("offset\tsum_cov\tn_covered\n")
    for j in range(2 * HALF + 1):
        out.write(f"{j - HALF}\t{sum_cov[j]}\t{n_cov[j]}\n")

if __name__ == "__main__":
    main()
