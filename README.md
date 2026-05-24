# cfDNA re-analysis — Snyder et al. 2016

Re-analysis of the cell-free DNA dataset from:

> Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. **Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin.** *Cell.* 2016 Jan 14;164(1-2):57–68. doi:10.1016/j.cell.2015.11.050. PMID: 26771485. GEO: [GSE71378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378).

## Study overview

The paper demonstrates that cfDNA fragmentation patterns in blood plasma reflect the genome-wide in vivo nucleosome occupancy of the cells that gave rise to those fragments, enabling inference of tissue-of-origin without relying on genotypic differences between tissues.

The dataset comprises **60 plasma cfDNA samples** across four groups:

- **Healthy (n=4):** BH01 is a pooled library from peripheral blood plasma contributed by an unknown number of healthy individuals (STEMCELL Technologies); IH01–IH03 are libraries from three individual healthy donors.
- **Autoimmune (n=8):** two donors with Crohn's disease (IA01–02), two with Ulcerative Colitis (IA03–04), and four with Systemic Lupus Erythematosus (IA05–08). Plasma was sourced from Conversant Bio or PlasmaLab International.
- **Cancer (n=48):** anonymous donors with clinical diagnoses of Stage IV solid tumours spanning over 20 cancer types. Forty-eight patients were enrolled, screened with light sequencing; four failed QC and were excluded (not deposited in GEO). Five of the remaining samples (IC15, IC17, IC20, IC35, IC37) were selected for deeper sequencing based on evidence of high tumour burden. Four additional cancer samples (IC49–IC52) appear in the GEO deposit, used in the combined CA01 nucleosome map.

Plasma cfDNA was purified with the QiaAMP Circulating Nucleic Acids kit (Qiagen) and quantified by Qubit fluorometer. Libraries were prepared with either a conventional **double-stranded protocol** (DSP; ThruPLEX-FD/DNA-seq, Rubicon Genomics) or a **single-stranded protocol** (SSP; adapted from Gansauge and Meyer 2013), the latter recovering short and damaged cfDNA fragments more efficiently. All libraries were sequenced on the Illumina HiSeq 2000. Reads were aligned to GRCh37 (1000 Genomes phase 2 reference) using BWA v0.7.10 (ALN algorithm), with PE reads further processed by BWA SAMPE.

## Files

### `snyder2016_metadata_GSE.csv`

The primary sample manifest for this re-analysis. Contains one row per sample (60 total) with all available metadata merged from three sources:

- **mmc1.pdf Table S4** — clinical data for cancer samples IC01–IC48: clinical diagnosis, stage, cfDNA yield (ng/ml), sex, QC outcome, and flags for samples selected for additional sequencing (`deep_seq`) or prepared from half the standard plasma volume (`half_vol_plasma`).
- **mmc1.pdf Tables S1 and S5** — sequencing statistics for all samples: library type (SSP/DSP), read configuration, fragment count, alignment rate, Q30 rate, coverage, estimated duplication rate, and proportion of fragments in the 35–80 bp (short/S-WPS) and 120–180 bp (long/L-WPS) size bins.
- **GEO series matrix (GSE71378)** — GEO sample accessions (`geo_accession`), BioSample accessions (`biosample_accession`), SRA experiment accessions (`sra_accession`), and disease/sex annotations for the full 60-sample set including non-cancer groups absent from Table S4.

An additional SRR-level accession (`srr_accession`) and file size information (`bases`, `bytes`) were joined from the SRA Run Table (see below). There is a 1:1 correspondence between SRX and SRR for all samples in this dataset.

Key columns:

| Column | Description |
|---|---|
| `sample_id` | Author-assigned sample label (e.g. BH01, IA03, IC17) |
| `geo_accession` | GEO sample accession (GSMxxxxxxx) |
| `biosample_accession` | NCBI BioSample accession (SAMNxxxxxxxx) |
| `sra_accession` | SRA experiment accession (SRXxxxxxxx) |
| `srr_accession` | SRA run accession (SRRxxxxxxx); use this for `prefetch`/`fasterq-dump` |
| `sample_group` | healthy / autoimmune / cancer |
| `disease` | Disease label from GEO (includes autoimmune subtype for IA samples) |
| `clinical_dx` | Clinical and histological diagnosis from Table S4 (cancer samples only) |
| `sex` | Patient sex |
| `cfdna_yield_ng_ml` | cfDNA yield from 1.0 ml plasma (ng/ml); Table S4 |
| `qc_fail` | TRUE for samples that failed QC and were excluded from analysis (not in GEO) |
| `deep_seq` | TRUE for the 5 samples selected for additional deep sequencing (§ in Table S4) |
| `half_vol_plasma` | TRUE for IC46, IC48 — only 0.5 ml plasma was available (\*\* in Table S4) |
| `library_type` | SSP (single-stranded) or DSP (double-stranded) |
| `read_config` | Sequencing read configuration; symmetric runs use `NxM` notation, asymmetric (R1 ≠ R2) use `N/M` notation (e.g. `43/42`) |
| `fragments_sequenced` | Total fragments sequenced after adapter trimming |
| `coverage` | Mean genome coverage |
| `pct_35_80bp` | % fragments in 35–80 bp range (short fraction, S-WPS input) |
| `pct_120_180bp` | % fragments in 120–180 bp range (long fraction, L-WPS input) |
| `bases` | Total bases in the SRA run |
| `bytes` | Compressed file size in the SRA run |

### `SraRunTable.csv`

Downloaded from the NCBI SRA Run Selector at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378. Contains one row per SRR run with NCBI-side metadata. The columns relevant to this project (`Run`, `Bases`, `Bytes`) have been merged into `snyder2016_metadata_GSE.csv`; this file is retained as the authoritative source for any additional SRA-level fields.

### `snyder2016_metadata_GSE_bam_subset.csv`

Three-sample subset of the main manifest used for alignment benchmarking: IH03 (SSP, 2×39, healthy), IA07 (DSP, 2×101, autoimmune), and IC25 (SSP, 43/42 asymmetric, cancer). These samples were chosen to cover both library protocols, symmetric and asymmetric read configurations, and a range of read lengths, and are used to compare the original Snyder alignments (downloaded BAMs) against re-alignment with modern tools.

## Scripts

### `download_fastq.py`

Downloads raw FASTQ files for all samples via `fasterq-dump` (SRA Toolkit). Reads `srr_accession` from an input CSV, skips samples whose FASTQs already exist, and reports per-sample timing and file sizes.

```
python3 download_fastq.py --input snyder2016_metadata_GSE.csv \
    --outdir data/fastq --tmpdir data/tmp --threads 6
```

All 60 samples were submitted to SRA as aligned BAMs; `fasterq-dump` converts them back to raw reads on the fly.

### `download_bam.py`

Downloads original aligned BAMs via `sam-dump | samtools view -bS`. Produces one `{SRR}.bam` per sample in the output directory. Intended for the benchmarking subset, but accepts any CSV with an `srr_accession` column.

```
python3 download_bam.py --input snyder2016_metadata_GSE_bam_subset.csv \
    --outdir data/bam --threads 6
```

Both scripts support background execution via `nohup ... > logfile.log 2>&1 &` and skip already-completed samples on re-run.

### `utils/extract_fastq_from_bam.py`

Extracts paired FASTQs from locally downloaded BAMs using `samtools collate | samtools fastq`. Preferred over `fasterq-dump` when BAMs are already on disk, as it avoids a redundant SRA download. Skips samples whose FASTQs already exist and skips BAMs that fail `samtools quickcheck` (i.e. are still downloading or truncated).

```
python3 utils/extract_fastq_from_bam.py \
    --input snyder2016_metadata_GSE_bam_subset.csv \
    --bamdir data/bam --outdir data/fastq --tmpdir data/tmp --threads 6
```

## Pipelines

### `cfdna-finale-snakemake/`

A Snakemake pipeline for end-motif and fragment-length analysis of cfDNA BAMs, built on [FinaleToolkit](https://epifluidlab.github.io/FinaleToolkit/). Implements the research questions in `End-motif-analysis.md`.

**Inputs:** original Snyder 2016 BAMs (`data/bam/`), an hg19 2-bit reference, and `samplesheet.csv`.

**Per-sample outputs** are produced for five fragment-length classes (`all`, `I`/`II`: sub- vs. mono-nucleosomal split, `III`/`IV`: Snyder S-WPS and L-WPS windows):

- 4-mer end-motif frequency tables and genome-wide MDS scores
- Per-100 kb interval end-motifs and MDS (BED)
- Full fragment-length histograms (1 bp bins) and distribution plots

Each rule runs in its own conda environment resolved via `--use-conda`; the driver environment (`environment.yml`) contains only Snakemake and pandas.

```bash
# on host, from cfdna-finale-snakemake/
bash launch.sh   # nohup; PID → snakemake.pid; log → snakemake.log
bash stop.sh
```

See `cfdna-finale-snakemake/README.md` for full configuration details.

### `ichorCNA/`

A Snakemake pipeline for tumor-fraction and copy-number estimation from low-coverage cfDNA BAMs using [ichorCNA](https://github.com/broadinstitute/ichorCNA). Applied to 56 samples (48 cancer IC02–IC52 + 8 autoimmune IA01–IA08); healthy normals are excluded and no panel of normals is used.

**Two rules per sample:**

1. **`read_counter`** — bins the BAM into 500 kb read-depth windows (WIG) with `readCounter` (HMMcopy suite), MAPQ ≥ 20.
2. **`ichorcna`** — runs the ichorCNA R script; GC-content and mappability WIGs are resolved from inside the conda-installed package, so no external reference files are needed.

**Key output:** `results/{sample}/{sample}.params.txt` — the `Tumor Fraction` field is the estimated ctDNA fraction. Additional outputs include copy-number segment files and genome-wide CN plots (PDF).

```bash
# on host, from ichorCNA/
bash launch.sh
bash stop.sh
```

See `ichorCNA/README.md` for parameter details and setup notes.
