# End-Motif Analysis of Snyder's cfDNA Data

## Motivation

End motifs — the 4-mer nucleotide sequence at each 5′ fragment end — reflect the identity of the nucleases responsible for cfDNA fragmentation rather than chromatin architecture per se. They are an orthogonal fragmentomic dimension relative to fragment size, WPS, and nucleosome footprints. See Jiang et al. 2020 for the foundational methodology and Snyder et al. 2016 for the short/long fragment framework this analysis builds on.

While fragment size and end motifs have been shown to co-vary (Jin et al. 2021), **no published study has applied end motif profiling within Snyder's specific S-WPS (35–120 bp) and L-WPS (120–220 bp) size classes, nor interpreted this stratification explicitly through the lens of differential nuclease activity in plasma cfDNA.** This is a concrete gap: short sub-nucleosomal fragments and long mono-nucleosomal fragments likely arise from different cell type and/or nuclease activities (DNASE1L3, DNASE1, caspase-activated DNase), and their end-motif landscapes should differ accordingly.

## Tool

**FinaleToolkit** (Li et al. 2025) is the practical standard. It takes BAM/CRAM as input and implements end motifs, MDS, interval-level end motifs, WPS, DELFI, and cleavage ratios in a single pip-installable Python package. The only comparable tools are Freefly (end motifs only, results essentially identical to FinaleToolkit: χ² = 0.001) and cfDNApro (R, OOM on large BAMs).

> Note: FinaleToolkit has **no** `bam-to-frags` / `frag.gz` conversion step. The native subcommands are `coverage, frag-length-bins, frag-length-intervals, cleavage-profile, wps, adjust-wps, delfi, delfi-gc-correct, end-motifs, interval-end-motifs, breakpoint-motifs, interval-breakpoint-motifs, mds, interval-mds, filter-file, agg-bw, gap-bed`. Size/quality filtering is done with `filter-file`, which emits a filtered **BAM**; all motif tools consume BAMs directly.

## Input Data

| Sample | Cohort | Library | Coverage | Reads |
| --- | --- | --- | --- | --- |
| IH02 | Single healthy individual | Single-stranded (ancient-DNA protocol) | 30× | 779 M |
| One of 44 cancer patients | Stage IV, various cancers | Single-stranded, light | 2.2× median | — |
| One of 5 cancer patients (selected) | SCLC, squamous lung, CRC, HCC, breast DCIS | Single-stranded | ~30× | — |

We should check Snyder 2016 / GSE71378 supplementary for reported tumor fraction and pick the highest-tumor-fraction samples. This gates the cancer-vs-normal questions and should be resolved before sample selection.

## Pipeline (validated Snakemake implementation)

Implemented and validated end-to-end on IH02 in `cfdna-finale-snakemake/`. The
commands below are the underlying FinaleToolkit invocations the Snakemake rules
issue; the workflow orchestrates them per sample and per length class. Key
parameters live in `config.yaml` (`mapq: 30`, `kmer: 4`,
`window_size_bp: 100000`, `length_classes`).

### 1. Reference preparation (contig-safe)

The Snyder BAMs are b37/1000G-style: bare contig names (`1..22, X, Y`), mito
`MT`, plus `hs37d5`/`NC_007605`/`GL000*.1` decoys. UCSC hg19 refs are
chr-prefixed, mito `chrM`, and carry alt-haplotype contigs (e.g.
`chr6_ssto_hap7`). FinaleToolkit builds intervals from the supplied 2-bit/
chrom.sizes and `fetch()`es each contig from the BAM — **any ref contig absent
from the BAM header crashes with `ValueError: invalid contig`.** The pipeline
therefore restricts refs to the intersection with the BAM `@SQ` header:

```bash
# authoritative contig list = BAM @SQ headers
samtools view -H sample.bam | awk '$1=="@SQ"' | sed 's/.*SN://; s/\t.*//' > bam_contigs.txt
# chrom.sizes restricted to BAM contigs
sed 's/^chr//' hg19.chrom.sizes | grep -Fwf bam_contigs.txt - > hg19_nochr.chrom.sizes
# 2bit subset to the same contigs (2bit ⊆ chrom.sizes ⊆ BAM @SQ)
cut -f1 hg19_nochr.chrom.sizes | sed 's/^/chr/' \
  | twoBitToFa -seqList=/dev/stdin hg19.2bit /dev/stdout \
  | sed 's/>chr/>/' | faToTwoBit /dev/stdin hg19_nochr.2bit
```

Net result for these BAMs: `1..22, X, Y`. MT, decoys, and alt contigs are
dropped — correct for cfDNA fragmentomics anyway.

### 2. Preprocess BAM

```bash
samtools sort -o sample.sorted.bam sample.bam
picard MarkDuplicates I=sample.sorted.bam O=sample.md.bam M=metrics.txt \
  CREATE_INDEX=false VALIDATION_STRINGENCY=LENIENT
samtools index sample.md.bam
```

### 3. Filter into size classes (`filter-file` → BAM)

`filter-file` auto-enforces: properly-paired, read1, not secondary/
supplementary, not QC-flagged — so MarkDuplicates-flagged duplicates are
dropped here. Size classing is done **at this step**, not inside `end-motifs`.

```bash
# all  (no length bounds — unstratified baseline)
finaletoolkit filter-file -q 30                  -w 4 -o sample.all.bam sample.md.bam
# I    < 120 (sub-nucleosomal, no lower bound)
finaletoolkit filter-file -q 30          -max 119 -w 4 -o sample.I.bam   sample.md.bam
# II   >= 120 (mono + di/tri-nucleosomal, open-ended tail)
finaletoolkit filter-file -q 30 -min 120          -w 4 -o sample.II.bam  sample.md.bam
# III  35-80 (Snyder strict S-WPS)
finaletoolkit filter-file -q 30 -min 35  -max 80  -w 4 -o sample.III.bam sample.md.bam
# IV   120-180 (Snyder strict L-WPS)
finaletoolkit filter-file -q 30 -min 120 -max 180 -w 4 -o sample.IV.bam  sample.md.bam
samtools index sample.{all,I,II,III,IV}.bam
```

`length_classes` (config): `all {–,–}`, `I {–, 119}`, `II {120, –}`,
`III {35, 80}`, `IV {120, 180}` — keyed directly to the analytical
questions. Classes overlap by design (I⊃III, II⊃IV); they answer different
questions, not a partition.

### 4. End motifs by size class

`end-motifs` requires the contig-safe 2-bit reference from step 1 as a
positional argument. `-v` enables a tqdm progress bar to the log.

```bash
for cls in all I II III IV; do
  finaletoolkit end-motifs -v -k 4 -q 30 -w 4 \
    -o sample.$cls.end_motifs.tsv  sample.$cls.bam  hg19_nochr.2bit
done
```

### 5. Motif Diversity Score

`mds` reads positionally and writes to **stdout**.

```bash
for cls in all I II III IV; do
  finaletoolkit mds sample.$cls.end_motifs.tsv > sample.$cls.mds.txt
done
```

MDS is the normalized Shannon entropy across all 256 4-mer frequencies. It
summarizes how evenly the 256 4-mers are distributed; lower MDS is associated
with DNASE1L3-dominated fragmentation (typical of healthy plasma) and higher
MDS with flatter distributions, but it can shift in either direction with
disease context (see Implementation notes).

### 6. Spatially resolved end motifs (per 100 kb window)

Li et al. 2025 showed interval-level MDS markedly outperforms genome-wide
summary MDS for early detection (AUC 0.44 → 0.84). Run per size class.

```bash
bedtools makewindows -g hg19_nochr.chrom.sizes -w 100000 > genome.100000bp.bed

for cls in all I II III IV; do
  finaletoolkit interval-end-motifs -v -k 4 -q 30 -w 4 \
    -o sample.$cls.interval_end_motifs.tsv \
    sample.$cls.bam  hg19_nochr.2bit  genome.100000bp.bed
  # interval-mds takes TWO positionals: file_path file_out
  finaletoolkit interval-mds sample.$cls.interval_end_motifs.tsv \
    sample.$cls.interval_mds.bed
done
```

### 7. Fragment-length distribution (QC / context)

```bash
finaletoolkit frag-length-bins -v -q 30 --bin-size 1 -min 1 -max 1000 \
  -o sample.frag_length_bins.tsv  sample.md.bam
```

`-min 1 -max 1000` is **required**: without length bounds a few discordant
read pairs push `max(TLEN)` to ~1e8, and FinaleToolkit's O(n_bins × n_keys)
binning loop then runs for days (see Implementation notes). The plot script
(`plot_frag_lengths.py`) reads the `min,max,count` schema (fragment length =
bin midpoint) and renders the full distribution plus 35–80 bp and 120–180 bp
sub-panels.

## Key Analytical Questions

Fragment length classes referenced below:
- I  < 120
- II ≥ 120
- III 35–80
- IV 120–180

> Classes overlap by design (I⊃III, II⊃IV): they target different questions, not a partition of fragment lengths.

1. Do 4-mers in I–IV differ between cancer and normal?
2. Do we see a difference in 2× vs 30× cancers (similar TF, same cancer type)?
3. Do end motifs differ between short and long fragments in cancer samples?

Additional questions:
- Do short and long fragments show distinct 4-mer end-motif profiles in healthy cfDNA?
- Is CCCA (a DNASE1L3-associated motif) enriched in long or short fragments, or distributed uniformly?
- Does MDS differ between length classes in healthy individuals, and does this shift in disease?
- Can size-stratified interval MDS improve upon genome-wide summary or size-pooled interval MDS?

## References

- Snyder 2016 — source data (GSE71378); established short/long fragment framework; single-stranded library protocol note
- Jiang et al. 2020 — foundational end-motif methodology; DNASE1L3/CCCA mechanism; MDS definition; size-stratified MDS precedent (on their own data, not Snyder's)
- FinaleToolkit — Li 2025 — implementation; interval-level end motif extension; benchmarking on BH01
