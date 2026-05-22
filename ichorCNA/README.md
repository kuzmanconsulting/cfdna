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
├── {sample}.wig              # read-depth bins from readCounter (step 1)
├── {sample}.params.txt       # estimated tumor fraction + ploidy
├── {sample}.seg              # copy-number segments
├── {sample}.cna.seg          # UCSC browser-compatible segments
├── {sample}.correctedDepth.txt
└── {sample}_*.pdf            # genome-wide CN plots
```

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

Uses the same Snakemake driver environment as the main pipeline
(`cfdna-finale-snakemake/environment.yml`). Must be run from the host, not
the container.

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
