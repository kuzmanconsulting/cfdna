# End-Motif Analysis of Snyder's cfDNA Data

## Motivation

End motifs — the 4-mer nucleotide sequence at each 5′ fragment end — reflect the identity of the nucleases responsible for cfDNA fragmentation rather than chromatin architecture per se. They are an orthogonal fragmentomic dimension relative to fragment size, WPS, and nucleosome footprints. See Jiang et al. 2020 for the foundational methodology and Snyder et al. 2016 for the short/long fragment framework this analysis builds on.

While fragment size and end motifs have been shown to co-vary (Jin et al. 2021), **no published study has applied end motif profiling within Snyder's specific S-WPS (35–120 bp) and L-WPS (120–220 bp) size classes, nor interpreted this stratification explicitly through the lens of differential nuclease activity in plasma cfDNA.** This is a concrete gap: short sub-nucleosomal fragments and long mono-nucleosomal fragments likely arise from different cell type and/or nuclease activities (DNASE1L3, DNASE1, caspase-activated DNase), and their end-motif landscapes should differ accordingly.

## Tool

**FinaleToolkit** (Li et al. 2025) is the practical standard. It takes BAM/CRAM or its native `frag.gz` format as input and implements end motifs, MDS, interval-level end motifs, WPS, DELFI, and cleavage ratios in a single pip-installable Python package. The only comparable tools are Freefly (end motifs only, results essentially identical to FinaleToolkit: χ² = 0.001) and cfDNApro (R, OOM on large BAMs).

## Input Data

| IH02                                    | Single healthy individual                  | Single-stranded (ancient DNA protocol) | 30×         | 779 M |
| --------------------------------------- | ------------------------------------------ | -------------------------------------- | ----------- | ----- |
| **One of** 44 cancer patients           | Stage IV, various cancers                  | Single-stranded, light                 | 2.2× median |       |
| **One of** 5 cancer patients (selected) | SCLC, squamous lung, CRC, HCC, breast DCIS | Single-stranded                        | ~30×        |       |

We should check if Snyder reports cancer fraction and use the highest cancer fraction  tumor samples.

## Proposed Pipeline

### 1. Preprocess BAM

```bash
# Sort and index (if not already done)
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam

# Mark duplicates
picard MarkDuplicates \
  I=sample.sorted.bam \
  O=sample.md.bam \
  M=metrics.txt
samtools index sample.md.bam
```

### 2. Convert to Fragment Format

```bash
finaletoolkit bam-to-frags \
  --input sample.md.bam \
  --output sample.frag.gz \
  --min-mapq 30 \
  --paired-end \
  --no-duplicates
```

Filters applied: MAPQ > 30, properly paired, duplicates excluded. The `frag.gz` format is more efficient for repeated analyses.

### 3. Compute End Motifs by Fragment Size Class

```bash
# Full pool
finaletoolkit end-motifs \
  --input sample.frag.gz \
  --output sample.end_motifs.all.tsv

# Short fragments (sub-nucleosomal; inclusive cutoff for motif count)
finaletoolkit end-motifs \
  --input sample.frag.gz \
  --min-length 35 --max-length 120 \
  --output sample.end_motifs.short.tsv

# Long fragments (mono-nucleosomal)
finaletoolkit end-motifs \
  --input sample.frag.gz \
  --min-length 120 --max-length 220 \
  --output sample.end_motifs.long.tsv
```

**On size cutoffs:** Snyder's original definitions (S-WPS: 35–80 bp; L-WPS: 120–180 bp) are strict and leave a coverage gap. For end-motif purposes, use broader ranges (short: <120 bp; long: 120–220 bp, analogous to DELFI's 100–150/151–220 split) to maximize fragment counts. Stricter cutoffs can be applied in sensitivity analyses.

Plot distribution of fragment lengths in full samples. Also plot fragment lengths distributions separately for 35-80bp and 120-180bp fragment lengths brackets.

### 4. Motif Diversity Score

```bash
finaletoolkit mds --input sample.end_motifs.all.tsv
finaletoolkit mds --input sample.end_motifs.short.tsv
finaletoolkit mds --input sample.end_motifs.long.tsv
```

MDS is the normalized Shannon entropy across all 256 4-mer frequencies. It summarizes nuclease diversity: low MDS = DNASE1L3-dominated (healthy); high MDS = flatter distribution (cancer, DNASE1L3 downregulation).

### 5. Spatially Resolved End Motifs (optional extension)

Li et al. 2025 demonstrated that **interval-level** MDS (per 100 kb window) dramatically outperforms genome-wide summary MDS for early cancer detection (AUC 0.44 → 0.84). Apply the same stratification here:

```bash
finaletoolkit interval-end-motifs \
  --input sample.frag.gz \
  --intervals 100kb_windows.bed \
  --output sample.interval_end_motifs.tsv

finaletoolkit interval-mds \
  --input sample.interval_end_motifs.tsv \
  --output sample.interval_mds.tsv
```

This should also be run separately per size class.

## Key Analytical Questions

Fragment length classes
- I  < 120
- II >= 120
- III 35-80
- IV 120-180

1. Do 4-mers in I-IV differ between cancer and normal?
2. Do we see difference in 2x vs 30x cancers? (similar TF and same cancer type)
3. Do we see difference in end motifs between short and long fragments in cancer samples?

Additional questions:
- Do short and long fragments show distinct 4-mer end-motif profiles in healthy cfDNA?
- Is CCCA (the DNASE1L3 signature motif) enriched in long or short fragments, or distributed uniformly?
- Does the MDS differ between length classes in healthy individuals, and does this difference shift in disease?
- Can size-stratified interval MDS improve upon the genome-wide summary or size-pooled interval MDS?

## Known issues

- **Size bin overlap at 120 bp.** Short is defined as 35–120 and long as 120–220, so a 120 bp fragment lands in both bins. 
long: over 120; short: 119 Finale toolkit should have an automated threshold
- **The `100kb_windows.bed` file just appears.** Step 5 calls it without saying where it comes from. FinaleToolkit doesn't ship it — we need to generate it (e.g., from a genome `chromsizes` file with `bedtools makewindows`).
- **MDS interpretation.** "Low MDS = DNASE1L3-dominated (healthy)" is correct directionally, but MDS can shift in both directions depending on disease context — it's really a measure of how evenly distributed the 256 4-mers are, not a direct DNASE1L3 proxy.

## References

- Snyder 2016 — source data (GSE71378); established short/long fragment framework; single-stranded library protocol note
- Jiang et al. 2020 — foundational end-motif methodology; DNASE1L3/CCCA mechanism; MDS definition; size-stratified MDS precedent (on their own data, not Snyder's)
- FinaleToolkit - Li 2025 — implementation; interval-level end motif extension; benchmarking on BH01
