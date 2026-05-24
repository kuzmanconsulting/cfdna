# ichorCNA

Tumor fraction and copy-number estimation from cfDNA BAMs using
[ichorCNA](https://github.com/broadinstitute/ichorCNA), applied to the
Snyder 2016 dataset.

## Inputs

- `config.yaml` — pipeline settings (see Parameters below)
- `samplesheet.csv` — one row per sample with `sample_id`, `srr_accession`,
  `sample_group`. BAMs are looked up as `{bam_dir}/{srr_accession}.bam`.

**Samples included:** 56 — 48 cancer (IC02–IC52) + 8 autoimmune (IA01–IA08).
Healthy normals (BH01, IH01–IH03) are excluded; no panel of normals is used.

## Outputs (per sample)

```
results/{sample}/
├── {sample}.params.txt            # tumor fraction, ploidy, all solutions
├── {sample}.seg.txt               # CNA segments — full table (use for analysis)
├── {sample}.seg                   # CNA segments — IGV/GISTIC-compatible
├── {sample}.cna.seg               # per-bin copy number and logR
├── {sample}.correctedDepth.txt    # GC/mappability-corrected log2 depth per bin
├── {sample}.wig                   # raw read counts per bin (readCounter output)
├── {sample}.RData                 # full R workspace (all HMM objects)
└── {sample}/                      # plots
    ├── {sample}_genomeWide.pdf         # genome-wide CNA profile — best solution
    ├── {sample}_genomeWide_all_sols.pdf# all initialisation solutions overlaid
    ├── {sample}_genomeWide_n*-p*.pdf   # per-initialisation genome-wide profiles
    ├── {sample}_CNA_chr*.pdf           # per-chromosome zoom-in (one per autosome)
    ├── {sample}_correct.pdf            # GC/mappability correction diagnostics
    ├── {sample}_bias.pdf               # depth vs. GC/mappability QC plots
    └── {sample}_tpdf.pdf               # HMM emission (t-distribution) fits
```

### File descriptions

**Data files**

| File | Description |
|------|-------------|
| `{sample}.params.txt` | **Start here.** Final tumor fraction and ploidy for the best-fit solution, plus GC-MAD (fit quality metric) and a table of all initialisation solutions with log-likelihoods. |
| `{sample}.seg.txt` | Segment-level CNA table: coordinates, bin count, median logR, integer copy-number call, subclonal status, and GC-corrected copy number. Most complete format — use for downstream analysis. Not IGV-compatible. |
| `{sample}.seg` | Same segments as `.seg.txt` in a simplified format compatible with IGV and GISTIC. |
| `{sample}.cna.seg` | Bin-level (500 kb resolution) copy-number calls and logR — one row per bin rather than per segment. Use when bin-level resolution is needed. |
| `{sample}.correctedDepth.txt` | GC- and mappability-corrected log2 read depth per 500 kb bin — the signal fed into the HMM, before segmentation. |
| `{sample}.wig` | Raw read counts per 500 kb bin in WIG format. Direct output of `readCounter`; input to ichorCNA. |
| `{sample}.RData` | Full R workspace saved after the run. Contains all intermediate HMM objects; needed for re-running with different parameters without recomputing WIGs. |

**Plots**

| File | Description |
|------|-------------|
| `{sample}_genomeWide.pdf` | **Main results plot.** Genome-wide copy-number profile for the best solution, annotated with tumor fraction and ploidy. |
| `{sample}_genomeWide_all_sols.pdf` | All initialisations overlaid — use to judge whether the chosen solution is clearly best or ambiguous. |
| `{sample}_genomeWide_n*-p*.pdf` | Genome-wide profile for each individual initialisation (normal fraction × ploidy combination). |
| `{sample}_CNA_chr*.pdf` | Chromosome-level zoom-in: bin logR with segment calls overlaid. One PDF per autosome. |
| `{sample}_correct.pdf` | Loess fit used for GC and mappability correction. Check if the genome-wide profile looks systematically biased. |
| `{sample}_bias.pdf` | Depth vs. GC content and mappability scatter, before and after correction. |
| `{sample}_tpdf.pdf` | Student's t mixture fits for each HMM emission state (neutral, gain, loss). |

The key result is `{sample}.params.txt`: the `Tumor Fraction` field is the
estimated ctDNA fraction for that sample.

## Pipeline

Two rules per sample:

1. **`read_counter`** — runs `readCounter` (HMMcopy suite) to bin the BAM into
   fixed-width read-depth windows (WIG format). 500 kb bins, MAPQ ≥ 20.

2. **`ichorcna`** — runs the ichorCNA R script. GC content and mappability WIG
   files are resolved at runtime from inside the conda env
   (`system.file("extdata", ..., package="ichorCNA")`), so no external
   reference files are needed.

## Parameters

### Alignment / binning

| Parameter | Value | Notes |
|---|---|---|
| `params.bin_size` | 500000 | readCounter window (bp). 500 kb is standard for ~2x cfDNA; finer bins require higher coverage. |
| `params.mapq` | 20 | MAPQ threshold for readCounter. ichorCNA convention is 20 (lower than the 30 used elsewhere in this project) to maximise read yield at low coverage. |
| `params.chrs` | `c(1:22)` | Autosomes only, no `chr` prefix — matches the GRCh37 / 1000G BAM headers. Must be consistent across readCounter and the ichorCNA R script. |

### ichorCNA HMM

| Parameter | Value | Notes |
|---|---|---|
| `params.normal_states` | `c(0.5,...,0.99)` | Grid of candidate normal fractions the HMM tries. ichorCNA picks the value maximising likelihood. The high end (0.95–0.99) covers the low tumor-fraction typical of cfDNA. |
| `params.ploidy` | `c(2,3)` | Candidate ploidies. The model is run at each ploidy independently; the best-fitting solution wins. Covers diploid and near-triploid tumors. |
| `params.max_cn` | 5 | Maximum callable copy-number state. Capped at 5 because depth noise at ~2x prevents reliable calling above this. |
| `params.include_homd` | FALSE | Whether to call homozygous deletions (CN=0). Off by default — HOMDs are rare and depth noise creates spurious zero-depth bins at shallow coverage. |
| `params.txne` | 0.9999 | HMM self-transition probability (probability of staying in the same CN state from one 500 kb bin to the next). Higher = longer, smoother segments. ichorCNA default. |
| `params.txn_strength` | 10000 | Scales the penalty on copy-number transitions. Works together with `txne` to control segmentation smoothness; leave both at defaults unless segments are visibly over- or under-fragmented. |

## Setup & running

Must be run from the host, not the container.

**Runtime dependencies are not managed by conda.** Although `use-conda: true`
is set in the Snakemake profile, neither the `ichorcna` nor `read_counter`
rules carry a `conda:` directive — they rely on tools already present on the
host system:

| Tool | Expected on host |
|------|-----------------|
| `Rscript` + `ichorCNA` R package | system R |
| `readCounter` | HMMcopy suite, in `$PATH` |
| `python3` + `matplotlib` | conda env `workflow/envs/plots.yaml` (auto-built by Snakemake) |

Install the ichorCNA R package once if not present:
```r
install.packages("BiocManager")
BiocManager::install("ichorCNA")   # or devtools::install_github("broadinstitute/ichorCNA")
```

```bash
micromamba activate snakemake-env
cd $HOST_WORKSPACE/ichorCNA
bash launch.sh   # nohup; PID in snakemake.pid; log in snakemake.log
bash stop.sh
```

### Verifying GC / mappability WIG names

After the conda envs are built (first run), confirm the in-package WIG names
match the `hg19_500kb` pattern assumed by the pipeline:

```bash
Rscript -e 'list.files(system.file("extdata", package="ichorCNA"), pattern="hg19")'
```

If the names differ, update the `_gc_wig()` and `_map_wig()` helper functions
in `workflow/rules/ichorcna.smk`.
