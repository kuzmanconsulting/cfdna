# Cell-of-origin cfDNA coverage analysis

Infers tissue-of-origin from the Snyder 2016 cfDNA by measuring read coverage at
**tissue-specific super-enhancers**. Active super-enhancers are nucleosome-depleted open
chromatin; when a tissue sheds cfDNA, fragments spanning its active regulatory regions are
preferentially cut, so coverage there drops. A **low cancer / matched-normal coverage ratio at
a tissue's super-enhancers therefore flags that tissue as a likely cell-of-origin.**

Super-enhancer coordinates come from **SEdb 3.0** (hg19), which matches the GRCh37/hg19 BAMs —
no liftOver. Coverage is profiled in a single **±1 kb window around each super-enhancer center**,
reduced to one scalar per sample × tissue, scaled by library size, and compared cancer-vs-normal.

## Route

| Step | Script | Runs on | Needs |
|---|---|---|---|
| Refs | `refs/download_refs_coo.sh` | host (network) | curl, awk |
| 0 | `00_build_sites.sh` | host | awk |
| 1 | `01_coverage.sh` | host | sambamba, samtools |
| 2 | `02_plot_coverage.R` | host | R + tidyverse |

`ROOT` defaults to `/mnt/scratch/DM/cfdna` (the host mount) in every script; override with
`ROOT=… bash …` if needed.

## Conda environment

`sambamba`/`samtools` for step 1 (the analysis container has neither):

```bash
micromamba create -f analysis/cell_of_origin/environment.yml   # genomics-tools
```

R / tidyverse for step 2 is handled by the existing R library, not this env.

## Config files

| File | Role |
|---|---|
| `coverage_samplesheet.csv` | the samples to profile: `sample_id, sample_group, library_type, coverage, tissues`. `tissues=all` measures the sample against every site set; a `;`-separated list restricts it. |
| `analysis/tumor_types.csv` | `disease, label, tissue` — the **single source of truth** mapping each Snyder diagnosis to its SEdb tissue (`Colorectal→colon`, `Lung→lung`, `Healthy→blood`, …). Drives both which reference tissues to download and each sample's own tissue. |

---

## Refs — `refs/download_refs_coo.sh`

Fetches the SEdb reference and the super-enhancer BEDs for exactly the tissues this panel needs.

### Process

1. **SEdb sample table** — downloads `Human_sample_information_sedb3.txt` (cols `Sample ID,
   Species, Data source, Biosample type, Tissue type, Biosample name`) to
   `refs/coo/sedb3_human_sample_information.txt` if absent.
2. **Panel tissues** — joins `coverage_samplesheet.csv` → `snyder2016_metadata_GSE.csv` (disease)
   → `tumor_types.csv` (tissue), giving the unique tissue set (here: blood, breast, colon, liver,
   lung, pancreas).
3. **Normal-lineage reference samples** — for each tissue, selects SEdb samples whose `Tissue
   type` matches the tissue **and** `Biosample type` ∈ {`Tissue`, `Primary cell`} — i.e. normal
   lineage, **cancer cell lines excluded**. The tissue→SEdb-vocabulary regex defaults to the
   tissue name as a substring; only `breast` (`breast|mammary`) and `liver` (`liver|hepat`) widen
   it. Result → `refs/coo/se_samples.tsv` (`sample_id, tissue`).
4. **Super-enhancer BEDs** — fetches each selected sample's whole-super-enhancer BED
   (`SE_hg19_bed/{id}_SE_hg19.bed`; cols `se_chr, se_start, se_end, se_id, …`) into
   `refs/coo/se_hg19/`. Cached — re-runs skip existing files; a non-`se_chr` header marks a bad
   download and is removed. We fetch SE, **not** `SE_ele`: coverage is profiled at the
   super-enhancer center, so only the SE coordinates are needed.

### Outputs

- `refs/coo/sedb3_human_sample_information.txt` — the SEdb sample table.
- `refs/coo/se_samples.tsv` — `sample_id, tissue` for the selected normal-lineage references.
- `refs/coo/se_hg19/{id}_SE_hg19.bed` — per-sample super-enhancer BEDs (cached).

### Run

```bash
bash refs/download_refs_coo.sh
```

---

## 0. Build sites — `00_build_sites.sh`

Turns each tissue's pooled super-enhancers into non-overlapping ±1 kb coverage windows that match
the BAMs. hg19-native, awk only — no liftOver, no bedtools.

### Inputs

| Input | What it is |
|---|---|
| `refs/coo/se_samples.tsv` | tissue → reference sample list (from refs step) |
| `refs/coo/se_hg19/*.bed` | per-sample super-enhancer BEDs (from refs step) |
| `refs/hg19.chrom.sizes` | bounds-checking windows against contig ends |
| params | `WINDOW=1000`, `MINGAP=2000` |

### Process — per tissue

1. **Pool centers** — for each reference sample of the tissue, emit one point per super-enhancer =
   its **midpoint** `int((se_start+se_end)/2)`.
2. **Standardize** — keep `chr1–22/X/Y`, **strip `chr`** → bare contigs matching the BAM `@SQ`
   names (b37), `sort -u`.
3. **Isolate** — greedily keep a center only if it's ≥ `MINGAP` (2 kb) from the last kept one on
   that chrom, so the ±1 kb windows never overlap. This also collapses near-duplicate
   super-enhancers pooled across samples.
4. **Expand + bound-check** — window = `[center−1000, center+1001)`, kept only if `start ≥ 0` and
   `end ≤ chrom size`.

### Outputs

- `sites/{tissue}.bed` — `chrom, start, end, center` (bare hg19 contigs, isolated 2001-bp windows).
- `sites/manifest.tsv` — `tissue, n_samples, n_superenhancers, n_sites` (sanity check).

### Run

```bash
micromamba run -n genomics-tools bash analysis/cell_of_origin/00_build_sites.sh \
    > analysis/cell_of_origin/build_sites.log 2>&1
```

---

## 1. Coverage — `01_coverage.sh`

For each sample × tissue, the mean cfDNA coverage in that tissue's ±1 kb windows, scaled by
library size. One scalar per pair. **Runs on host** — sambamba is not in the analysis container.

### Inputs

| Input | What it is |
|---|---|
| `sites/{tissue}.bed` | per-tissue ±1 kb windows from step 0 |
| `coverage_samplesheet.csv` | which samples to run, and which tissues per sample |
| `cfdna-finale-snakemake/work/{sample}/{sample}.md.bam` (+ `.bai`) | preprocessed BAMs: coordinate-sorted, MarkDuplicates-flagged, indexed |
| env / params | `SAMBAMBA`/`SAMTOOLS` (default on `PATH`), `THREADS` (default 8); filter = `MAPQ≥30, not duplicate/secondary/qcfail/unmapped` |

### Process — per sample × tissue

1. **Resolve BAM** — `…/work/{sample}/{sample}.md.bam`; skip the sample if missing. Dedup flags
   and the index must already exist (that's why this uses the pipeline's `work/` BAMs).
2. **Library size** — `samtools idxstats` → sum of mapped reads = `total_reads`.
3. **Mean coverage** — `sambamba depth region` over the tissue's windows with the MAPQ/dedup
   filter; mean of the `meanCoverage` column across windows = `mean_cov`.
4. **Scale** — `scaled_cov = mean_cov / total_reads × 1e6`. The constant is cosmetic — it cancels
   in the cancer/normal ratio, which is what corrects for differing library depths.
5. **Atomic, idempotent write** — `> .tmp` then `mv`; existing outputs are skipped.

The `tissues` column gates which site sets a sample is measured against (`all` = every tissue).

### Outputs

- `coverage/{sample}.{tissue}.cov.tsv` — one row: `sample, tissue, n_windows, mean_cov,
  total_reads, scaled_cov`. Consumed by `02_plot_coverage.R`.

### Run

```bash
micromamba run -n genomics-tools bash analysis/cell_of_origin/01_coverage.sh \
    > analysis/cell_of_origin/coverage.log 2>&1
```

---

## 2. Plot — `02_plot_coverage.R`

Computes the cancer / matched-normal super-enhancer coverage ratio and plots it. Run from repo
root (relative paths).

### Process

1. **Load** all `coverage/*.cov.tsv` (one row per sample × tissue).
2. **Annotate** from `snyder2016_metadata_GSE.csv` + `tumor_types.csv`: `sample_group`, the
   sample's `own_tissue`, and its matched-control **stratum**:
   - `DSP → IH01`, `SSP & coverage ≥ 10 → IH02`, `SSP & coverage < 10 → IH03`.
3. **Ratio** — for each cancer sample × tissue, `ratio = scaled_cov / scaled_cov[matched normal,
   same tissue]`; `log2ratio`; flag `is_own` where `tissue == own_tissue`.

### Outputs

- `coo_coverage_ratio.tsv` — the full per-(cancer sample, tissue) ratio table.
- `coo_coverage_heatmap.png` — log2(cancer/normal) heatmap, sample × tissue, own-tissue cells
  outlined in black. Negative (blue) at a sample's own tissue = cell-of-origin signal.
- `coo_own_tissue.png` — barplot of the log2 ratio at each cancer sample's own tissue.

### Run

```bash
Rscript analysis/cell_of_origin/02_plot_coverage.R
```
