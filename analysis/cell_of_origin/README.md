# Cell-of-origin cfDNA coverage analysis

## Conda environment

```bash
micromamba create -f analysis/cell_of_origin/environment.yml
```

## 0. Build sites

A per-tissue loop that turns two hg38 databases into hg19 coverage-profiling windows matching your BAMs.

### Inputs

| Input | What it is |
|---|---|
| `refs/coo/gtrd_hg38_tf_centers.bed.gz` | ~51 M GTRD TF-binding-site centers (hg38): `chrom, pos, pos+1, TF` |
| `refs/coo/sedb3_human_sample_information.txt` | SEdb sample table: `Sample ID, …, Biosample type, Tissue type, …` |
| SEdb per-sample SE BEDs | downloaded on demand from `…/SE_hg38_bed/{id}_SE_hg38.bed` (hg38 super-enhancer regions) |
| `refs/coo/liftOver` + `hg38ToHg19.over.chain.gz` | UCSC liftOver to convert hg38→hg19 |
| `refs/hg19.chrom.sizes` | for bounds-checking windows against contig ends |
| `tissue_map.tsv` | `tissue → SEdb Tissue-type regex` (colon, lung, liver, breast, pancreas, blood) |
| params | `WINDOW=1000`, `MINGAP=2000`, `MAXSAMPLES=60`, `MAXSITES=20000` |

### Process — for each tissue (lines 38–96)

1. **Select samples** (43–45): from the SEdb table, take Sample IDs whose `Tissue type` matches the tissue regex **and** `Biosample type` is `Tissue`/`Primary cell` — i.e. normal-lineage SEs, *not* cancer cell lines. Cap at `MAXSAMPLES`.

2. **Download + pool SEs** (49–58): fetch each sample's SE BED (cached in `se_dl/`), verify the `se_chr` header, drop it, keep cols 1–3 (`chrom/start/end`), append to a raw pool.

3. **SE union** (61–63): filter to standard chroms (`chr1–22/X/Y`), `sort`, `bedtools merge` → the tissue's non-overlapping super-enhancer regions (hg38). *This is the cell-type filter.*

4. **Intersect GTRD ∩ SE union** (66–69): `bedtools intersect -u` keeps every GTRD center landing inside any tissue SE → the TF binding sites active in that tissue's regulome (hg38). *GTRD supplies precise point anchors; SE membership makes them tissue-specific.*

5. **liftOver hg38→hg19** (72–73): lift centers to hg19 so they match the BAMs; unmapped centers are dropped.

6. **Finalize in hg19** (77–91):
   - **strip `chr`** + keep `1–22/X/Y` → bare contigs matching the BAM `@SQ` names;
   - **sort** by position;
   - **isolate** (greedy, line 79): keep a center only if it's ≥ `MINGAP` (2 kb) from the last kept one — so the ±1 kb windows never overlap. This is what lets `aggregate_offsets.py` map each base to exactly one center;
   - **cap** (83–84, 88): if > `MAXSITES` remain, take every *step*-th to bound coverage runtime;
   - **expand + bound-check** (89–90): window = `[center−1000, center+1001)`, kept only if `start≥0` and `end ≤ chrom size`.

### Outputs

- **`sites/{tissue}.bed`** — `chrom, start, end, center` (bare hg19 contigs, isolated ±1 kb windows). Consumed by `01_coverage.sh`; the `center` column drives the aggregator's offset math.
- **`sites/manifest.tsv`** — per-tissue counts: `n_samples, n_se_union, n_centers_hg38, n_sites_final` (your sanity check).
- **`se_tmp/`** — intermediates (raw pools, SE unions, hg38/hg19 centers, isolated lists) + cached SE downloads in `se_tmp/se_dl/` (so re-runs skip re-downloading).

**Net transformation:** *"all TF binding sites genome-wide" (hg38) → "tissue-active TF binding sites, spaced ≥2 kb apart, as hg19 ±1 kb windows."* The spacing and hg19/bare-contig conversion exist to make the downstream coverage aggregation correct against the BAMs.

### Run

```bash
nohup micromamba run -n genomics-tools \
    bash analysis/cell_of_origin/00_build_sites.sh \
    > analysis/cell_of_origin/build_sites.log 2>&1
```

## 1. Coverage

For each sample × tissue, computes per-base cfDNA coverage in the tissue's ±1 kb windows and collapses it to one −1000…+1000 mean-coverage profile.

### Inputs

| Input | What it is |
|---|---|
| `sites/{tissue}.bed` | per-tissue isolated ±1 kb windows from step 0 (`chrom, start, end, center`) |
| `cfdna-finale-snakemake/work/{sample}/{sample}.md.bam` (+ `.bai`) | preprocessed BAMs: coordinate-sorted, MarkDuplicates-flagged, indexed |
| `cfdna-finale-snakemake/samplesheet.csv` | the 18-sample list (`sample_id, srr, …`) |
| `aggregate_offsets.py` | companion script that collapses per-base depth to an offset profile |
| env / params | `SAMBAMBA` (default `sambamba`), `THREADS` (default 8); filter = `MAPQ≥30, not duplicate/secondary/qcfail/unmapped`, all fragment sizes |

### Process — for each sample × tissue (`01_coverage.sh` 33–45, `aggregate_offsets.py`)

1. **Discover tissues** (29–31): list `sites/*.bed` → tissue names.

2. **Resolve BAM** (33–36): the preprocessed `…/work/{sample}/{sample}.md.bam`; skip the sample if it's missing. Dedup flags and the index must already exist — that's why this uses the pipeline's `work/` BAMs, not raw `data/bam/`.

3. **Skip-if-done** (37–40): `coverage/{sample}.{tissue}.offset_cov.tsv` is skipped if present, so the step is idempotent and resumable.

4. **`sambamba depth base`** (41): per-base depth over that tissue's windows (`-L`), applying the MAPQ/dedup `-F` filter; emits `chrom pos cov …` for each covered base.

5. **Collapse to offsets** (`aggregate_offsets.py`): load the sites BED → per-chrom sorted center list + total `n_sites` (py 23–32); stream each `(chrom, pos, cov)`, binary-search the single window containing `pos` (isolated ⇒ ≤1 match), compute `offset = pos − center`, and accumulate `sum_cov[offset]` + `n_covered[offset]` (py 39–68); write the profile (py 70–74).

6. **Atomic write** (42–43): `> .tmp` then `mv`, so a killed run never leaves a partial profile.

### Outputs

- **`coverage/{sample}.{tissue}.offset_cov.tsv`** — a `# n_sites=… sample=… tissue=…` comment line, then `offset, sum_cov, n_covered` for offsets −1000…+1000. **Mean coverage at an offset = `sum_cov / n_sites`** — the denominator is the window count, *not* `n_covered`, because zero-coverage bases are simply absent from sambamba's stream. Consumed by `02_plot_coverage.R`.

The large depth spread across samples (~2×–148×) is **not** corrected here — that's handled by flank-normalization in step 2; this step stores raw summed coverage.

### Run

```bash
micromamba run -n genomics-tools \
    bash analysis/cell_of_origin/01_coverage.sh \
    > analysis/cell_of_origin/coverage.log 2>&1
```
