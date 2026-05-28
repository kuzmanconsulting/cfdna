#!/usr/bin/env Rscript
# Load samples from cfdna-finale-snakemake/samplesheet.csv, attach disease
# (from snyder2016_metadata_GSE.csv) and disease_label (from
# analysis/tumor_types.csv).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggridges)
  library(patchwork)
  library(yaml)
})

config <- yaml::read_yaml("analysis/config.yaml")

dir.create("analysis/all_spls",  showWarnings = FALSE, recursive = TRUE)
dir.create("analysis/ssp_subset", showWarnings = FALSE, recursive = TRUE)

log_con <- file("analysis/analysis.log", open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) })

# Shared palettes
group_colors <- c(healthy = "#523FCC", cancer = "#FF004C")
lib_colors   <- c(SSP = "#2166AC", DSP = "#B2182B")

samples      <- read_csv("cfdna-finale-snakemake/samplesheet.csv", show_col_types = FALSE)
meta         <- read_csv("snyder2016_metadata_GSE.csv",            show_col_types = FALSE)
tumor_types  <- read_csv("analysis/tumor_types.csv",               show_col_types = FALSE)

samples <- samples |>
  left_join(meta |> select(sample_id, disease), by = "sample_id") |>
  left_join(tumor_types, by = "disease") |>
  rename(disease_label = label) |>
  mutate(library_type = factor(library_type, levels = c("SSP", "DSP")),
         sample_group = factor(sample_group, levels = c("healthy", "cancer")))

missing <- samples |> filter(is.na(disease_label))
if (nrow(missing))
  stop("No disease_label for: ", paste(unique(missing$disease), collapse = "; "))

print(samples, n = Inf)

# Coverage vs tumor fraction overview
p_cov_tf <- ggplot(samples, aes(tumor_fraction, coverage,
                                shape = library_type, color = sample_group)) +
  geom_vline(xintercept = 0.1, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_point(size = 2.2) +
  geom_text_repel(aes(label = if_else(sample_group == "healthy",
                                       sample_id,
                                       paste(sample_id, disease_label))),
                  size = 2.6, color = "grey40", alpha = 0.7,
                  min.segment.length = 0, segment.size = 0.2,
                  segment.color = "grey70", max.overlaps = Inf) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Tumor fraction (ichorCNA)", y = "Coverage (log10)",
       shape = "Library", color = "Group") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("analysis/all_spls/coverage_vs_tf.png", p_cov_tf,
       width = 5, height = 4, dpi = 300)
message("\nWrote analysis/all_spls/coverage_vs_tf.png")

# ---------------------------------------------------------------------------
# Fragment length ridgeplot (35-500 bp)
# ---------------------------------------------------------------------------
cache_dir  <- "analysis/.cache"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cache_file <- file.path(cache_dir, "frag_lengths_35_500.tsv")

if (!file.exists(cache_file)) {
  message("Building frag-length cache: ", cache_file)
  frag <- map(samples$sample_id, \(s) {
    f <- sprintf("cfdna-finale-snakemake/results/%s/frag_lengths/%s.frag_length_bins.tsv", s, s)
    read_tsv(f, show_col_types = FALSE) |>
      transmute(sample_id = s, length = min, count = count) |>
      filter(length >= 35, length <= 500)
  }) |> list_rbind()
  write_tsv(frag, cache_file)
} else {
  message("Using cached frag-length file: ", cache_file)
  frag <- read_tsv(cache_file, show_col_types = FALSE)
}

frag <- frag |>
  left_join(samples |> select(sample_id, sample_group, library_type,
                              disease_label, coverage),
            by = "sample_id") |>
  mutate(cov_tier = if_else(coverage >= 10, "ge10", "lt10"))

sample_order <- samples |>
  mutate(group_rank = if_else(sample_group == "healthy", 0L, 1L)) |>
  arrange(group_rank, desc(coverage)) |>
  mutate(y_label = paste0(sample_id, " (", disease_label, ")"))

frag <- frag |>
  left_join(sample_order |> select(sample_id, y_label), by = "sample_id") |>
  group_by(cov_tier) |>
  mutate(height = count / max(count)) |>   # scale within each tier
  ungroup()

make_ridge <- function(tier, title) {
  d <- frag |> filter(cov_tier == tier)
  lvls <- rev(sample_order |>
                filter(sample_id %in% d$sample_id) |>
                pull(y_label))
  d <- d |> mutate(y_label = factor(y_label, levels = lvls))
  ggplot(d, aes(x = length, y = y_label, height = height,
                fill = library_type, color = library_type)) +
    annotate("rect", xmin = 35,  xmax = 80,  ymin = -Inf, ymax = Inf,
             fill = "#FDD0A2", alpha = 0.35) +
    annotate("rect", xmin = 120, xmax = 180, ymin = -Inf, ymax = Inf,
             fill = "#C7E9C0", alpha = 0.35) +
    geom_vline(xintercept = 120, linetype = "dashed",
               color = "grey30", linewidth = 0.3) +
    geom_ridgeline(scale = 0.9, alpha = 0.35, linewidth = 0.4) +
    scale_fill_manual(values = lib_colors, limits = names(lib_colors), drop = FALSE) +
    scale_color_manual(values = lib_colors, limits = names(lib_colors), drop = FALSE) +
    scale_x_continuous(limits = c(35, 500), breaks = seq(50, 500, 50),
                       expand = c(0, 0)) +
    labs(title = title, x = "Fragment length (bp)", y = NULL,
         fill = "Library", color = "Library") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
}

n_top <- sum(samples$coverage >= 10)
n_bot <- sum(samples$coverage <  10)

p_top <- make_ridge("ge10", "Coverage ≥ 10×")
p_bot <- make_ridge("lt10", "Coverage < 10×") +
  guides(fill = "none", color = "none")

p_frag <- (p_top / p_bot) +
  plot_layout(guides = "collect", heights = c(n_top, n_bot)) &
  theme(legend.position = "right")

ggsave("analysis/all_spls/frag_length_ridges.png", p_frag,
       width = 8, height = 0.35 * nrow(samples) + 2, dpi = 300)
message("Wrote analysis/all_spls/frag_length_ridges.png")

# ---------------------------------------------------------------------------
# Median 4-mer end-motif barplots (per class × library/group combo)
# ---------------------------------------------------------------------------
motif_cache <- file.path(cache_dir, "end_motifs_I_IV.tsv")
classes <- c("I", "II", "III", "IV")

if (!file.exists(motif_cache)) {
  message("Building end-motif cache: ", motif_cache)
  motifs <- tidyr::expand_grid(sample_id = samples$sample_id, class = classes) |>
    pmap(\(sample_id, class) {
      f <- sprintf("cfdna-finale-snakemake/results/%s/end_motifs/%s.%s.tsv",
                   sample_id, sample_id, class)
      read_tsv(f, col_names = c("kmer", "freq"), col_types = "cd") |>
        mutate(sample_id = sample_id, class = class)
    }) |> list_rbind()
  write_tsv(motifs, motif_cache)
} else {
  message("Using cached end-motif file: ", motif_cache)
  motifs <- read_tsv(motif_cache, show_col_types = FALSE)
}

motifs <- motifs |>
  left_join(samples |> select(sample_id, sample_group, library_type, coverage),
            by = "sample_id") |>
  mutate(class      = factor(class, levels = classes),
         first_nuc  = substr(kmer, 1, 1),
         kmer       = factor(kmer, levels = sort(unique(kmer))))

nuc_colors <- c(A = "#4393C3", T = "#FFC107", C = "#D6604D", G = "#4DAC26")

make_motif_bar <- function(d, facet_x) {
  med <- d |>
    group_by(class, kmer, first_nuc, across(all_of(facet_x))) |>
    summarise(median_freq = median(freq), .groups = "drop")

  ggplot(med, aes(kmer, median_freq, fill = first_nuc)) +
    geom_col(width = 1, linewidth = 0) +
    scale_fill_manual(values = nuc_colors, name = "First nt") +
    scale_x_discrete(breaks = NULL) +
    facet_grid(reformulate(facet_x, "class")) +
    labs(x = "4-mer motif (256, sorted alphabetically)",
         y = "Median frequency") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey92"))
}

# Top: >=10x, faceted by library × group (4 cols)
top_d <- motifs |>
  filter(coverage >= 10) |>
  mutate(lib_group = factor(paste(library_type, sample_group, sep = "_"),
                            levels = c("SSP_healthy", "SSP_cancer",
                                       "DSP_healthy", "DSP_cancer")))

p_top_mot <- make_motif_bar(top_d, "lib_group") +
  ggtitle("Median 4-mer end-motif frequencies, coverage ≥ 10×")

ggsave("analysis/all_spls/motif_freq_ge10.png", p_top_mot,
       width = 12, height = 7, dpi = 300)
message("Wrote analysis/all_spls/motif_freq_ge10.png")

# Bottom: <10x, faceted by group only (SSP only — no DSP samples below 10×)
bot_d <- motifs |> filter(coverage < 10)
stopifnot(all(bot_d$library_type == "SSP"))

p_bot_mot <- make_motif_bar(bot_d, "sample_group") +
  ggtitle("Median 4-mer end-motif frequencies, coverage < 10× (SSP only)")

ggsave("analysis/all_spls/motif_freq_lt10.png", p_bot_mot,
       width = 7, height = 7, dpi = 300)
message("Wrote analysis/all_spls/motif_freq_lt10.png")

# ---------------------------------------------------------------------------
# MDS (normalized Shannon entropy) strip plot
# ---------------------------------------------------------------------------
mds_cache <- file.path(cache_dir, "mds_I_IV_all.tsv")

if (!file.exists(mds_cache)) {
  message("Building MDS cache: ", mds_cache)
  read_mds <- function(s, cls) {
    f <- sprintf("cfdna-finale-snakemake/results/%s/end_motifs/%s.%s.mds.txt",
                 s, s, cls)
    if (!file.exists(f)) { warning("missing: ", f); return(NA_real_) }
    as.numeric(readLines(f, n = 1))
  }
  mds_tbl <- tidyr::expand_grid(sample_id = samples$sample_id,
                                class = c(classes, "all")) |>
    mutate(mds = map2_dbl(sample_id, class, read_mds))
  write_tsv(mds_tbl, mds_cache)
} else {
  message("Using cached MDS file: ", mds_cache)
  mds_tbl <- read_tsv(mds_cache, show_col_types = FALSE)
}

mds_tbl <- mds_tbl |>
  left_join(samples |> select(sample_id, sample_group, library_type,
                              disease_label, coverage),
            by = "sample_id") |>
  mutate(label = if_else(sample_group == "healthy",
                         sample_id, paste(sample_id, disease_label)))

mds_class <- mds_tbl |>
  filter(class %in% classes) |>
  mutate(class = factor(class, levels = classes))

dodge_w <- 0.5
g_lvls  <- levels(mds_class$sample_group)
n_g     <- length(g_lvls)
offsets <- setNames(-dodge_w/2 + (seq_len(n_g) - 0.5) * dodge_w / n_g, g_lvls)

mds_class <- mds_class |>
  mutate(x_num = as.integer(class) + offsets[as.character(sample_group)])

p_mds <- ggplot(mds_class, aes(x_num, mds,
                               color = sample_group, shape = library_type)) +
  geom_point(aes(alpha = coverage), size = 2, stroke = 0.8) +
  scale_alpha_continuous(range = c(0.25, 1), trans = "log10") +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = c(SSP = 16, DSP = 17), guide = "none") +
  scale_x_continuous(breaks = seq_along(classes), labels = classes) +
  scale_y_continuous(limits = c(0.90, 1)) +
  facet_wrap(~ library_type, nrow = 1) +
  labs(x = "Length class", y = "Motif diversity score",
       color = "Group", alpha = "Coverage") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

ggsave("analysis/all_spls/mds_strip.png", p_mds,
       width = 8, height = 4, dpi = 300)
message("Wrote analysis/all_spls/mds_strip.png")

# Per-sample plots
if (isTRUE(config$skip_per_spl_4mer_freq)) {
  message("Skipping per-sample 4-mer plots (config$skip_per_spl_4mer_freq = TRUE)")
} else {
per_sample_dir <- "analysis/all_spls/4mer_freqs"
dir.create(per_sample_dir, showWarnings = FALSE, recursive = TRUE)

for (s in samples$sample_id) {
  d <- motifs |> filter(sample_id == s)
  dl <- samples$disease_label[samples$sample_id == s]
  lt <- samples$library_type[samples$sample_id == s]
  cv <- samples$coverage[samples$sample_id == s]
  p <- ggplot(d, aes(kmer, freq, fill = first_nuc)) +
    geom_col(width = 1, linewidth = 0) +
    scale_fill_manual(values = nuc_colors, name = "First nt") +
    scale_x_discrete(breaks = NULL) +
    facet_grid(class ~ .) +
    labs(title = sprintf("%s (%s), %s, %.1fX", s, dl, lt, cv),
         x = "4-mer motif (256, sorted alphabetically)",
         y = "Frequency") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey92"))
  ggsave(file.path(per_sample_dir, sprintf("%s_%s.png", s, lt)), p,
         width = 8, height = 7, dpi = 200)
}
message("Wrote per-sample 4-mer plots to ", per_sample_dir)
}

# ---------------------------------------------------------------------------
# Per-motif log2 fold-change: cancer (median) vs appropriate healthy control
#   a) DSP         → cancer IC49/50/51/52 vs IH01
#   b) SSP ≥10×    → cancer IC15/17/20/35/37 vs IH02
#   c) SSP <10×    → cancer IC10/23/28/32/33/46 vs IH03
# ---------------------------------------------------------------------------
message("\n--- Per-motif log2 fold-change ---")

# Classify every sample into one of three strata
motifs_strat <- motifs |>
  mutate(stratum = case_when(
    library_type == "DSP"                  ~ "DSP",
    library_type == "SSP" & coverage >= 10 ~ "SSP_ge10",
    library_type == "SSP" & coverage <  10 ~ "SSP_lt10"
  ))

strata_order  <- c("DSP", "SSP_ge10", "SSP_lt10")
stratum_ctrl  <- c(DSP = "IH01", SSP_ge10 = "IH02", SSP_lt10 = "IH03")
stratum_label <- c(DSP = "DSP", SSP_ge10 = "SSP ≥10×", SSP_lt10 = "SSP <10×")

# Report stratum membership
for (st in strata_order) {
  ids <- motifs_strat |> filter(stratum == st, sample_group == "cancer") |>
    pull(sample_id) |> unique() |> sort()
  message(sprintf("  %s  ctrl=%s  cancer: %s",
                  st, stratum_ctrl[st], paste(ids, collapse = ", ")))
}

# Compute log2FC per stratum
lfc_all <- map(strata_order, function(st) {
  d <- motifs_strat |> filter(stratum == st)

  ctrl_freq <- d |>
    filter(sample_group == "healthy") |>
    group_by(class, kmer) |>
    summarise(ctrl_freq = mean(freq), .groups = "drop")   # 1 ctrl sample; mean == identity

  cancer_freq <- d |>
    filter(sample_group == "cancer") |>
    group_by(class, kmer, first_nuc) |>
    summarise(cancer_median = median(freq), .groups = "drop")

  cancer_freq |>
    left_join(ctrl_freq, by = c("class", "kmer")) |>
    mutate(
      log2fc  = log2(cancer_median / ctrl_freq),
      stratum = st
    )
}) |>
  list_rbind() |>
  mutate(
    stratum = factor(stratum, levels = strata_order),
    class   = factor(class,   levels = classes)
  )

# Warn about any non-finite values
bad <- lfc_all |> filter(!is.finite(log2fc))
if (nrow(bad)) {
  warning(sprintf("%d non-finite log2FC values (zero freq in ctrl or cancer); dropped from plot",
                  nrow(bad)))
  lfc_all <- lfc_all |> filter(is.finite(log2fc))
}

make_lfc_bar <- function(st) {
  d       <- lfc_all  |> filter(stratum == st)
  ctrl_id <- stratum_ctrl[st]
  lbl     <- stratum_label[st]
  n_cancer <- motifs_strat |>
    filter(stratum == st, sample_group == "cancer") |>
    pull(sample_id) |> n_distinct()

  ggplot(d, aes(kmer, log2fc, fill = first_nuc)) +
    geom_col(width = 1, linewidth = 0) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey30") +
    scale_fill_manual(values = nuc_colors, name = "First nt") +
    scale_x_discrete(breaks = NULL) +
    facet_grid(class ~ ., scales = "free_y") +
    labs(
      title = sprintf("%s — cancer median (n=%d) vs %s (healthy control)",
                      lbl, n_cancer, ctrl_id),
      x = "4-mer motif (256, sorted alphabetically)",
      y = "log₂ FC (cancer / healthy)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid       = element_blank(),
      axis.ticks.x     = element_blank(),
      strip.background = element_rect(fill = "grey92")
    )
}

out_files <- c(DSP = "lfc_dsp.png", SSP_ge10 = "lfc_ssp_ge10.png", SSP_lt10 = "lfc_ssp_lt10.png")
for (st in strata_order) {
  p <- make_lfc_bar(st)
  ggsave(file.path("analysis/all_spls", out_files[st]), p, width = 8, height = 7, dpi = 300)
  message("Wrote analysis/all_spls/", out_files[st])
}

# ---------------------------------------------------------------------------
# PCA — all 18 samples, length classes I–IV concatenated as features
#
# Feature matrix: 18 samples × 1024 columns (256 kmers × 4 classes).
# Concatenation encodes fragment length structurally: short-class (I, III)
# and long-class (II, IV) motif vectors occupy distinct columns so that
# length-specific signals pull samples along different PC axes.
# center=TRUE, scale.=FALSE — features are already commensurate frequencies.
# ---------------------------------------------------------------------------
message("\n--- PCA: all samples, length classes I–IV concatenated ---")

# Wide feature matrix: one row per sample, one column per "kmer_class"
pca_wide <- motifs |>
  mutate(feature = paste0(as.character(kmer), "_", as.character(class))) |>
  select(sample_id, feature, freq) |>
  pivot_wider(names_from = feature, values_from = freq)

mat_pca <- pca_wide |> select(-sample_id) |> as.matrix()
rownames(mat_pca) <- pca_wide$sample_id
stopifnot(!anyNA(mat_pca))

pca_res <- prcomp(mat_pca, center = TRUE, scale. = FALSE)
pct_var  <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

message(sprintf("  PC1 = %.1f%%  PC2 = %.1f%%  PC3 = %.1f%%",
                pct_var[1], pct_var[2], pct_var[3]))

# Sample scores
scores_pca <- as_tibble(pca_res$x[, 1:2], rownames = "sample_id") |>
  left_join(samples |> select(sample_id, sample_group, library_type, coverage),
            by = "sample_id")

# Top loadings by vector magnitude in PC1–PC2 plane, scaled to sample score range
top_n_load <- 12L
load_scale <- max(abs(pca_res$x[, 1:2])) / max(abs(pca_res$rotation[, 1:2]))

top_idx <- order(
  sqrt(pca_res$rotation[, 1]^2 + pca_res$rotation[, 2]^2),
  decreasing = TRUE
)[seq_len(top_n_load)]

loadings_pca <- as_tibble(pca_res$rotation[top_idx, 1:2], rownames = "feature") |>
  mutate(
    across(c(PC1, PC2), \(x) x * load_scale),
    kmer  = substr(feature, 1, 4),              # always exactly 4 chars
    class = substr(feature, 6, nchar(feature))  # drop "KMER_"
  )

p_pca_all <- ggplot() +
  # Loading arrows (grey; class annotated in label)
  geom_segment(
    data      = loadings_pca,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow     = arrow(length = unit(0.12, "cm"), type = "closed"),
    linewidth = 0.35, color = "grey55", alpha = 0.85
  ) +
  geom_text_repel(
    data               = loadings_pca,
    aes(PC1, PC2, label = paste0(kmer, " (", class, ")")),
    size               = 2.1, color = "grey30",
    min.segment.length = 0, segment.size = 0.2, segment.color = "grey70",
    max.overlaps       = Inf, box.padding = 0.2
  ) +
  # Sample points
  geom_point(
    data = scores_pca,
    aes(PC1, PC2, color = sample_group, shape = library_type, alpha = coverage),
    size = 3, stroke = 0.7
  ) +
  geom_text_repel(
    data               = scores_pca,
    aes(PC1, PC2, label = sample_id, color = sample_group),
    size               = 2.6, fontface = "bold", show.legend = FALSE,
    min.segment.length = 0, segment.size = 0.2, segment.color = "grey70",
    max.overlaps       = Inf, box.padding = 0.35
  ) +
  scale_color_manual(values = group_colors, name = "Group") +
  scale_shape_manual(values = c(SSP = 16, DSP = 17), name = "Library") +
  scale_alpha_continuous(range = c(0.25, 1), trans = "log10", name = "Coverage") +
  labs(
    title = sprintf("End-motif PCA — all %d samples, classes I–IV concatenated (1024 features)",
                    nrow(scores_pca)),
    x = sprintf("PC1 (%.1f%%)", pct_var[1]),
    y = sprintf("PC2 (%.1f%%)", pct_var[2])
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("analysis/all_spls/pca_all_samples.png", p_pca_all, width = 7, height = 6, dpi = 300)
message("Wrote analysis/all_spls/pca_all_samples.png")

# ---------------------------------------------------------------------------
# PCA — per library type: SSP and DSP side by side (patchwork)
#
# Each PCA is fitted independently on its library subset so that PC axes
# reflect within-library variance only. The dominant between-library split
# (DSP vs SSP strand-selection effect) is removed by design, letting the
# cancer vs healthy signal occupy the principal axes.
# ---------------------------------------------------------------------------
message("\n--- PCA: per-library (SSP | DSP) ---")

top_n_load_lib <- 10L

# Shared coverage scale — identical limits + breaks → patchwork merges the legends
cov_limits <- c(1, 200)
cov_breaks <- c(2, 10, 30, 100)

make_lib_pca <- function(lib, show_tf_legend = TRUE) {
  d_wide <- motifs |>
    filter(as.character(library_type) == lib) |>
    mutate(feature = paste0(as.character(kmer), "_", as.character(class))) |>
    select(sample_id, feature, freq) |>
    pivot_wider(names_from = feature, values_from = freq)

  mat <- d_wide |> select(-sample_id) |> as.matrix()
  rownames(mat) <- d_wide$sample_id

  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  message(sprintf("  %s  PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%  n=%d",
                  lib, pct[1], pct[2], pct[3], nrow(mat)))

  scores <- as_tibble(pca$x[, 1:2], rownames = "sample_id") |>
    left_join(samples |> select(sample_id, sample_group, library_type,
                                coverage, tumor_fraction),
              by = "sample_id")

  load_scale <- max(abs(pca$x[, 1:2])) / max(abs(pca$rotation[, 1:2]))
  top_idx <- order(
    sqrt(pca$rotation[, 1]^2 + pca$rotation[, 2]^2),
    decreasing = TRUE
  )[seq_len(top_n_load_lib)]

  loadings <- as_tibble(pca$rotation[top_idx, 1:2], rownames = "feature") |>
    mutate(
      across(c(PC1, PC2), \(x) x * load_scale),
      kmer  = substr(feature, 1, 4),
      class = substr(feature, 6, nchar(feature))
    )

  ggplot() +
    geom_segment(
      data      = loadings,
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow     = arrow(length = unit(0.12, "cm"), type = "closed"),
      linewidth = 0.35, color = "grey55", alpha = 0.85
    ) +
    geom_text_repel(
      data               = loadings,
      aes(PC1, PC2, label = paste0(kmer, " (", class, ")")),
      size               = 2.1, color = "grey30",
      min.segment.length = 0, segment.size = 0.2, segment.color = "grey70",
      max.overlaps       = Inf, box.padding = 0.2
    ) +
    geom_point(
      data = scores,
      aes(PC1, PC2, color = sample_group, shape = library_type,
          alpha = coverage, size = tumor_fraction),
      stroke = 0.7
    ) +
    geom_text_repel(
      data               = scores,
      aes(PC1, PC2, label = sample_id, color = sample_group),
      size               = 2.6, fontface = "bold", show.legend = FALSE,
      min.segment.length = 0, segment.size = 0.2, segment.color = "grey70",
      max.overlaps       = Inf, box.padding = 0.35
    ) +
    scale_color_manual(values = group_colors, name = "Group") +
    scale_shape_manual(values = c(SSP = 16, DSP = 17), guide = "none") +
    scale_alpha_continuous(range = c(0.25, 1), trans = "log10", name = "Coverage",
                           limits = cov_limits, breaks = cov_breaks) +
    scale_size_continuous(range = c(1.5, 7), name = "Tumor fraction",
                          labels = scales::label_percent(accuracy = 1),
                          guide  = if (show_tf_legend) "legend" else "none") +
    labs(
      title = sprintf("%s  (n = %d)", lib, nrow(scores)),
      x     = sprintf("PC1 (%.1f%%)", pct[1]),
      y     = sprintf("PC2 (%.1f%%)", pct[2])
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

p_ssp <- make_lib_pca("SSP")
p_dsp <- make_lib_pca("DSP", show_tf_legend = FALSE)

p_pca_libs <- (p_ssp | p_dsp) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "End-motif PCA by library — classes I–IV concatenated (1024 features)"
  ) &
  theme(legend.position = "right")

ggsave("analysis/all_spls/pca_by_library.png", p_pca_libs, width = 14, height = 6, dpi = 300)
message("Wrote analysis/all_spls/pca_by_library.png")

# ---------------------------------------------------------------------------
# SSP subset: IC17, IC37, IC20, IC35 vs IH02
#   Criteria: SSP library, ~30× coverage, ~0.2 tumor fraction
#   Analyses:
#     A) Median 4-mer log2FC across the 4 cancer samples (all + I–IV)
#     B) Per-sample log2FC — full 256-kmer barplots (all + I–IV)
#     C) Concordance report + dot plot of fully-concordant motifs
# ---------------------------------------------------------------------------
message("\n--- SSP subset: IC17/IC37/IC20/IC35 vs IH02 ---")

ssp_sub_cancer  <- c("IC17", "IC37", "IC20", "IC35")
ssp_sub_ctrl    <- "IH02"
ssp_sub_ids     <- c(ssp_sub_ctrl, ssp_sub_cancer)
ssp_all_classes <- c("all", "I", "II", "III", "IV")

# ---- 1. Cache ----------------------------------------------------------
ssp_sub_cache <- file.path(cache_dir, "end_motifs_ssp_subset_all5.tsv")

if (!file.exists(ssp_sub_cache)) {
  message("Building SSP subset motif cache: ", ssp_sub_cache)
  ssp_sub_motifs <- tidyr::expand_grid(sample_id = ssp_sub_ids,
                                       class     = ssp_all_classes) |>
    pmap(\(sample_id, class) {
      f <- sprintf("cfdna-finale-snakemake/results/%s/end_motifs/%s.%s.tsv",
                   sample_id, sample_id, class)
      read_tsv(f, col_names = c("kmer", "freq"), col_types = "cd") |>
        mutate(sample_id = sample_id, class = class)
    }) |> list_rbind()
  write_tsv(ssp_sub_motifs, ssp_sub_cache)
} else {
  message("Using cached SSP subset motifs: ", ssp_sub_cache)
  ssp_sub_motifs <- read_tsv(ssp_sub_cache, show_col_types = FALSE)
}

ssp_sub_motifs <- ssp_sub_motifs |>
  mutate(
    class     = factor(class, levels = ssp_all_classes),
    first_nuc = substr(kmer, 1, 1),
    kmer      = factor(kmer, levels = sort(unique(kmer)))
  )

# ---- 2. Per-sample log2FC (cancer vs IH02) ----------------------------
ctrl_motifs <- ssp_sub_motifs |>
  filter(sample_id == ssp_sub_ctrl) |>
  select(class, kmer, ctrl_freq = freq)

persample_lfc_sub <- ssp_sub_motifs |>
  filter(sample_id %in% ssp_sub_cancer) |>
  left_join(ctrl_motifs, by = c("class", "kmer")) |>
  mutate(log2fc = log2(freq / ctrl_freq))

bad_sub <- persample_lfc_sub |> filter(!is.finite(log2fc))
if (nrow(bad_sub))
  warning(sprintf("%d non-finite log2FC values in SSP subset; dropped", nrow(bad_sub)))
persample_lfc_sub <- persample_lfc_sub |> filter(is.finite(log2fc))

# ---- A. Median log2FC barplot ------------------------------------------
median_lfc_sub <- persample_lfc_sub |>
  group_by(class, kmer, first_nuc) |>
  summarise(log2fc = median(log2fc), .groups = "drop")

p_ssp_median <- ggplot(median_lfc_sub, aes(kmer, log2fc, fill = first_nuc)) +
  geom_col(width = 1, linewidth = 0) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey30") +
  scale_fill_manual(values = nuc_colors, name = "First nt") +
  scale_x_discrete(breaks = NULL) +
  facet_grid(class ~ ., scales = "free_y") +
  labs(
    title = sprintf("SSP subset: median cancer log₂FC vs %s  (n=%d: %s)",
                    ssp_sub_ctrl, length(ssp_sub_cancer),
                    paste(ssp_sub_cancer, collapse = ", ")),
    x     = "4-mer motif (256, sorted alphabetically)",
    y     = "Median log₂ FC (cancer / IH02)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid       = element_blank(),
        axis.ticks.x     = element_blank(),
        strip.background = element_rect(fill = "grey92"))

ggsave("analysis/ssp_subset/lfc_ssp_subset_median.png", p_ssp_median,
       width = 8, height = 10, dpi = 300)
message("Wrote analysis/ssp_subset/lfc_ssp_subset_median.png")

# ---- B. Per-sample log2FC grid (4 samples × 5 classes) ----------------
# Build informative row-strip labels
sub_meta <- samples |>
  filter(sample_id %in% ssp_sub_cancer) |>
  mutate(
    strip_label = sprintf("%s\n%s\n%.0f×, TF≈%.0f%%",
                          sample_id, disease_label,
                          coverage, tumor_fraction * 100)
  ) |>
  arrange(match(sample_id, ssp_sub_cancer))   # keep consistent order

strip_levels <- sub_meta$strip_label

persample_lfc_sub2 <- persample_lfc_sub |>
  left_join(sub_meta |> select(sample_id, strip_label), by = "sample_id") |>
  mutate(strip_label = factor(strip_label, levels = strip_levels))

p_ssp_persample <- ggplot(persample_lfc_sub2,
                          aes(kmer, log2fc, fill = first_nuc)) +
  geom_col(width = 1, linewidth = 0) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey30") +
  scale_fill_manual(values = nuc_colors, name = "First nt") +
  scale_x_discrete(breaks = NULL) +
  facet_grid(strip_label ~ class, scales = "free_y") +
  labs(
    x = "4-mer motif (256, sorted alphabetically)",
    y = "log₂ FC (vs IH02)"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid        = element_blank(),
        axis.ticks.x      = element_blank(),
        strip.background  = element_rect(fill = "grey92"),
        strip.text.y      = element_text(size = 7.5))

ggsave("analysis/ssp_subset/lfc_ssp_subset_per_sample.png", p_ssp_persample,
       width = 12, height = 7, dpi = 300)
message("Wrote analysis/ssp_subset/lfc_ssp_subset_per_sample.png")

# ---- C. Concordance ----------------------------------------------------
concordance_sub <- persample_lfc_sub |>
  group_by(class, kmer, first_nuc) |>
  summarise(
    n_up       = sum(log2fc > 0),
    n_down     = sum(log2fc < 0),
    n_total    = n(),
    median_lfc = median(log2fc),
    sd_lfc     = sd(log2fc),
    .groups    = "drop"
  ) |>
  mutate(
    n_concord     = pmax(n_up, n_down),
    fully_concord = n_concord == length(ssp_sub_cancer),
    direction     = if_else(n_up >= n_down, "up", "down")
  )

message("\n  Concordant motifs (all 4 samples agree on direction):")
for (cls in levels(concordance_sub$class)) {
  hits <- concordance_sub |> filter(class == cls, fully_concord)
  n_fc <- nrow(hits)
  top_up <- hits |> filter(direction == "up")  |> arrange(desc(median_lfc)) |> head(8)
  top_dn <- hits |> filter(direction == "down") |> arrange(median_lfc)       |> head(8)
  message(sprintf("  class %-4s: %3d/256 fully concordant", cls, n_fc))
  if (nrow(top_up)) message("    Up:   ",
    paste(sprintf("%s(%.3f)", top_up$kmer, top_up$median_lfc), collapse = "  "))
  if (nrow(top_dn)) message("    Down: ",
    paste(sprintf("%s(%.3f)", top_dn$kmer, top_dn$median_lfc), collapse = "  "))
}

# ---- C-csv. Concordant motifs wide CSV ---------------------------------
concord_col_order <- c(
  paste0("fq_", c("all", "I", "II", "III", "IV"), "_lfc_dn"),
  paste0("fq_", c("all", "I", "II", "III", "IV"), "_lfc_up")
)

concord_csv <- concordance_sub |>
  filter(fully_concord) |>
  mutate(col_name = paste0("fq_", as.character(class), "_lfc_",
                           if_else(direction == "up", "up", "dn"))) |>
  select(kmer, col_name, median_lfc) |>
  pivot_wider(names_from = col_name, values_from = median_lfc)

# Ensure every expected column is present even when no motif is concordant
# in a particular class/direction
for (col in concord_col_order)
  if (!col %in% names(concord_csv)) concord_csv[[col]] <- NA_real_

concord_csv <- concord_csv |>
  rename(motif = kmer) |>
  select(motif, all_of(concord_col_order))

write_csv(concord_csv, "analysis/ssp_subset/ssp_subset_concordant_motifs.csv", na = "")
message(sprintf("Wrote analysis/ssp_subset/ssp_subset_concordant_motifs.csv  (%d motifs)",
                nrow(concord_csv)))

if (nrow(concord_csv) == 0) {
  message("  No fully concordant motifs found — skipping concordant dot plot")
} else {

# Concordant motif dot plot — top 15 up + top 15 down per class
top_concord <- concordance_sub |>
  filter(fully_concord) |>
  group_by(class, direction) |>
  slice_max(abs(median_lfc), n = 15, with_ties = FALSE) |>
  ungroup()

# Per-class ordered factor for y-axis (ordered by median_lfc within each class)
# Use "class::kmer" prefix so global factor levels respect per-class ordering
kmer_levels <- top_concord |>
  arrange(class, median_lfc) |>
  mutate(ck = paste0(as.character(class), "::", as.character(kmer))) |>
  pull(ck)

concord_pts <- persample_lfc_sub |>
  semi_join(top_concord, by = c("class", "kmer")) |>
  mutate(ck = factor(paste0(as.character(class), "::", as.character(kmer)),
                     levels = kmer_levels))

concord_med <- top_concord |>
  mutate(ck = factor(paste0(as.character(class), "::", as.character(kmer)),
                     levels = kmer_levels))

sample_colors_sub <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),   # 4 Set1-style hues
  ssp_sub_cancer
)
sample_labels_sub <- setNames(
  paste0(sub_meta$sample_id, " (", sub_meta$disease_label, ")"),
  sub_meta$sample_id
)

p_concord <- ggplot() +
  # zero line
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  # spine from 0 to median
  geom_segment(data = concord_med,
               aes(x = 0, xend = median_lfc, y = ck, yend = ck),
               linewidth = 0.5, color = "grey60") +
  # individual sample points
  geom_point(data  = concord_pts,
             aes(x = log2fc, y = ck, color = sample_id),
             size  = 2, alpha = 0.85) +
  # median marker
  geom_point(data  = concord_med,
             aes(x = median_lfc, y = ck),
             shape = 21, size = 3, fill = "grey25", color = "white", stroke = 0.5) +
  scale_color_manual(values = sample_colors_sub, labels = sample_labels_sub,
                     name = "Sample") +
  scale_y_discrete(labels = \(x) sub(".*::", "", x)) +
  facet_wrap(~ class, nrow = 1, scales = "free_y") +
  labs(
    title = sprintf(
      "Concordant 4-mer shifts vs IH02 — top 15 up/down per class (4/4 samples agree)"),
    x     = "log₂ FC (vs IH02)",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey92"),
        axis.text.y        = element_text(family = "mono", size = 8))

ggsave("analysis/ssp_subset/lfc_ssp_subset_concordant.png", p_concord,
       width = 14, height = 7, dpi = 300)
message("Wrote analysis/ssp_subset/lfc_ssp_subset_concordant.png")

} # end if (nrow(concord_csv) > 0)

# ---------------------------------------------------------------------------
# Class III vs IV concordant motif comparison
#   III = 35–80 bp (sub-nucleosomal), IV = 120–180 bp (mono-nucleosomal)
#   Question: which motifs shift concordantly in one size class but not the
#   other, and do any flip direction between short and long fragments?
# ---------------------------------------------------------------------------
message("\n--- Class III vs IV concordant motif comparison ---")

iii_iv <- concordance_sub |>
  filter(class %in% c("III", "IV")) |>
  select(class, kmer, first_nuc, median_lfc, fully_concord, direction) |>
  pivot_wider(
    names_from  = class,
    values_from = c(median_lfc, fully_concord, direction),
    names_sep   = "_"
  ) |>
  mutate(
    concord_III = replace_na(fully_concord_III, FALSE),
    concord_IV  = replace_na(fully_concord_IV,  FALSE),
    same_dir    = direction_III == direction_IV,
    membership  = case_when(
      concord_III & concord_IV & !same_dir ~ "opposite direction",
      concord_III & concord_IV             ~ "III ∩ IV (same dir)",
      concord_III                          ~ "III only",
      concord_IV                           ~ "IV only",
      TRUE                                 ~ "neither"
    ) |> factor(levels = c("opposite direction", "III ∩ IV (same dir)",
                           "III only", "IV only", "neither"))
  )

# Console report
n_iii     <- sum(iii_iv$concord_III)
n_iv      <- sum(iii_iv$concord_IV)
n_shared  <- sum(iii_iv$concord_III & iii_iv$concord_IV)
n_opp     <- sum(iii_iv$membership == "opposite direction")

message(sprintf("  III concordant : %d/256", n_iii))
message(sprintf("  IV  concordant : %d/256", n_iv))
message(sprintf("  Shared (III∩IV): %d  (%d same direction, %d opposite)",
                n_shared, n_shared - n_opp, n_opp))
message(sprintf("  III only       : %d", sum(iii_iv$membership == "III only")))
message(sprintf("  IV  only       : %d", sum(iii_iv$membership == "IV only")))

fmt_motifs <- function(d, lfc_col) {
  d |> arrange(desc(abs(.data[[lfc_col]]))) |>
    mutate(s = sprintf("%s(%.3f)", kmer, .data[[lfc_col]])) |>
    pull(s) |> paste(collapse = "  ")
}

for (grp in c("III only", "IV only", "III ∩ IV (same dir)", "opposite direction")) {
  d <- iii_iv |> filter(membership == grp)
  if (nrow(d) == 0) next
  lfc_col <- if (grp == "IV only") "median_lfc_IV" else "median_lfc_III"
  message(sprintf("\n  %s (%d motifs):", grp, nrow(d)))
  for (dir in c("up", "down")) {
    dd <- d |> filter(if (grp == "IV only") direction_IV == dir else direction_III == dir)
    if (nrow(dd)) message(sprintf("    %s: %s", dir, fmt_motifs(dd, lfc_col)))
  }
}

# Scatter: all 256 motifs, median_lfc_III vs median_lfc_IV
memb_colors <- c(
  "opposite direction"     = "#B2182B",
  "III ∩ IV (same dir)"   = "#907396ff",
  "III only"               = "#2166AC",
  "IV only"                = "#E08214",
  "neither"                = "grey82"
)
memb_sizes <- c(
  "opposite direction"    = 3,
  "III ∩ IV (same dir)"  = 2.5,
  "III only"              = 2,
  "IV only"               = 2,
  "neither"               = 1
)
memb_alpha <- c(
  "opposite direction"    = 1,
  "III ∩ IV (same dir)"  = 1,
  "III only"              = 0.9,
  "IV only"               = 0.9,
  "neither"               = 0.35
)

label_data_iii_iv <- iii_iv |> filter(membership != "neither")

p_iii_iv <- ggplot(iii_iv, aes(median_lfc_III, median_lfc_IV, color = membership)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, color = "grey60") +
  geom_point(aes(size = membership, alpha = membership)) +
  geom_text_repel(
    data               = label_data_iii_iv,
    aes(label         = as.character(kmer)),
    size               = 2.6, family = "mono",
    min.segment.length = 0, segment.size = 0.2, segment.color = "grey60",
    max.overlaps       = Inf, box.padding = 0.25
  ) +
  scale_x_continuous(limits = c(-0.7, 0.7)) +
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  scale_color_manual(values = memb_colors, name = NULL) +
  scale_size_manual(values  = memb_sizes,  guide = "none") +
  scale_alpha_manual(values = memb_alpha,  guide = "none") +
  labs(
    title = "Class III vs IV: concordant cancer motif shifts (SSP subset vs IH02)",
    x     = "Median log₂ FC — class III (35–80 bp)",
    y     = "Median log₂ FC — class IV (120–180 bp)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right")

ggsave("analysis/ssp_subset/lfc_ssp_subset_III_vs_IV.png", p_iii_iv,
       width = 8, height = 7, dpi = 300)
message("\nWrote analysis/ssp_subset/lfc_ssp_subset_III_vs_IV.png")

# ---- Grouped barplot: motifs discordant between class III and IV ----------
# Filter to motifs concordant in exactly one of III/IV (or opposite direction)
diff_motifs <- iii_iv |>
  filter(membership %in% c("III only", "IV only", "opposite direction"))

# Order within each membership by the primary class log2FC
kmer_order_bar <- c(
  diff_motifs |> filter(membership == "opposite direction") |>
    arrange(median_lfc_III) |> pull(kmer) |> as.character(),
  diff_motifs |> filter(membership == "III only") |>
    arrange(median_lfc_III) |> pull(kmer) |> as.character(),
  diff_motifs |> filter(membership == "IV only") |>
    arrange(median_lfc_IV)  |> pull(kmer) |> as.character()
)

# Pull median + SD for III and IV from concordance_sub; join membership
diff_bar <- concordance_sub |>
  filter(as.character(class) %in% c("III", "IV"),
         as.character(kmer)  %in% kmer_order_bar) |>
  left_join(diff_motifs |> select(kmer, membership), by = "kmer") |>
  mutate(
    kmer  = factor(as.character(kmer),  levels = kmer_order_bar),
    class = factor(as.character(class), levels = c("III", "IV"))
  )

class_colors_bar <- c(III = "#2166AC", IV = "#E08214")

p_db <- ggplot(diff_bar, aes(x = kmer, y = median_lfc, fill = class)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey40") +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, linewidth = 0) +
  geom_errorbar(
    aes(ymin = median_lfc - sd_lfc, ymax = median_lfc + sd_lfc),
    position  = position_dodge(width = 0.75),
    width     = 0.35, linewidth = 0.4, color = "grey30"
  ) +
  scale_fill_manual(values = class_colors_bar, name = "Class") +
  facet_grid(. ~ membership, scales = "free_x", space = "free_x") +
  labs(
    title   = "Class-discordant concordant shifts: III only vs IV only (SSP subset vs IH02)",
    x       = NULL,
    y       = "Median log₂ FC (vs IH02) ± SD",
    caption = "Concordance defined as all 4 cancer samples (IC17/IC37/IC20/IC35) agreeing on direction vs IH02"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                    family = "mono", size = 7),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey92"),
    legend.position    = "right",
    plot.caption       = element_text(size = 7.5, color = "grey45", hjust = 0,
                                      margin = margin(t = 6))
  )

p_db_width <- nrow(diff_motifs) * 0.22 + 2.5

ggsave("analysis/ssp_subset/lfc_ssp_subset_III_vs_IV_barplot.png", p_db,
       width = 10, height = 5, dpi = 300)
message("Wrote analysis/ssp_subset/lfc_ssp_subset_III_vs_IV_barplot.png")

# ---- Paired III − IV difference per motif ---------------------------------
# For each motif × cancer sample compute log2FC_III − log2FC_IV directly.
# Avoids the binary concordance threshold; measures class-specificity as a
# continuous paired quantity.

top_n_diff <- 20L

diff_paired <- persample_lfc_sub |>
  filter(as.character(class) %in% c("III", "IV")) |>
  select(sample_id, kmer, first_nuc, class, log2fc) |>
  pivot_wider(names_from = class, values_from = log2fc, names_prefix = "lfc_") |>
  mutate(diff_III_IV = lfc_III - lfc_IV)

diff_summary <- diff_paired |>
  group_by(kmer, first_nuc) |>
  summarise(median_diff = median(diff_III_IV),
            sd_diff     = sd(diff_III_IV),
            .groups = "drop")

top_diff <- bind_rows(
  diff_summary |> slice_max(median_diff, n = top_n_diff, with_ties = FALSE),
  diff_summary |> slice_min(median_diff, n = top_n_diff, with_ties = FALSE)
) |>
  distinct() |>
  arrange(median_diff) |>
  mutate(kmer = factor(as.character(kmer), levels = as.character(kmer)))

diff_pts <- diff_paired |>
  filter(as.character(kmer) %in% levels(top_diff$kmer)) |>
  mutate(kmer = factor(as.character(kmer), levels = levels(top_diff$kmer)))

p_diff_paired <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_segment(data = top_diff,
               aes(x = 0, xend = median_diff, y = kmer, yend = kmer),
               linewidth = 0.5, color = "grey65") +
  geom_point(data  = diff_pts,
             aes(x = diff_III_IV, y = kmer, color = sample_id),
             size = 2.2, alpha = 0.85) +
  geom_point(data  = top_diff,
             aes(x = median_diff, y = kmer, color = "Median"),
             shape = 19, size = 3.2) +
  scale_color_manual(values = c(sample_colors_sub, Median = "grey20"),
                     labels = c(sample_labels_sub, Median = "Median"),
                     name = "Sample") +
  labs(
    title = sprintf("Top %d motifs by class III − IV log₂FC (SSP subset vs IH02)", top_n_diff),
    x     = "log₂FC(III) − log₂FC(IV)  [per sample]",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(family = "mono", size = 9),
    legend.position    = "right"
  )

ggsave("analysis/ssp_subset/lfc_ssp_subset_III_vs_IV_diff.png", p_diff_paired,
       width = 7, height = top_n_diff * 0.28 + 2, dpi = 300)
message("Wrote analysis/ssp_subset/lfc_ssp_subset_III_vs_IV_diff.png")

# ---------------------------------------------------------------------------
# Focal motif log2FC: 6 motifs of interest across length classes I–IV
#   CCCA, CCAG, CCTG, TAAA, AAAA, TTTT
#   SSP subset cancer samples (IC17/IC37/IC20/IC35) vs IH02
# ---------------------------------------------------------------------------
message("\n--- Focal motif log2FC: selected motifs across classes I–IV ---")

focal_motifs   <- c("CCCA", "CCAG", "CCTG", "TAAA", "AAAA", "TTTT")
focal_classes  <- c("I", "II", "III", "IV")

focal_pts <- persample_lfc_sub |>
  filter(as.character(kmer) %in% focal_motifs,
         as.character(class) %in% focal_classes) |>
  mutate(
    kmer  = factor(as.character(kmer),  levels = focal_motifs),
    class = factor(as.character(class), levels = focal_classes)
  )

focal_med <- focal_pts |>
  group_by(kmer, class) |>
  summarise(median_lfc = median(log2fc), .groups = "drop")

# Colour each cancer sample distinctly (reuse palette from concordance section)
focal_sample_colors <- sample_colors_sub   # already defined above

p_focal <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey50") +
  # per-sample points + connecting lines across classes
  geom_line(data  = focal_pts,
            aes(x = class, y = log2fc, group = sample_id, color = sample_id),
            linewidth = 0.55, alpha = 0.6) +
  geom_point(data = focal_pts,
             aes(x = class, y = log2fc, color = sample_id),
             size = 2.2, alpha = 0.85) +
  # median line + marker
  geom_line(data  = focal_med,
            aes(x = class, y = median_lfc, group = kmer),
            linewidth = 0.9, color = "grey20", alpha = 0.9) +
  geom_point(data  = focal_med,
             aes(x = class, y = median_lfc),
             shape = 21, size = 3.5, fill = "grey20",
             color = "white", stroke = 0.6) +
  scale_color_manual(values = focal_sample_colors,
                     labels = sample_labels_sub,
                     name   = "Sample") +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  facet_wrap(~ kmer, nrow = 2, scales = "free_y") +
  labs(
    x = "Fragment length class",
    y = "log₂ FC (vs IH02)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background   = element_rect(fill = "grey92"),
    strip.text         = element_text(family = "mono", face = "bold", size = 10),
    legend.position    = "right"
  )

ggsave("analysis/ssp_subset/lfc_ssp_subset_focal_motifs.png", p_focal,
       width = 9, height = 5.5, dpi = 300)
message("Wrote analysis/ssp_subset/lfc_ssp_subset_focal_motifs.png")

# ---------------------------------------------------------------------------
# SSP subset MDS across all 5 fragment-length classes
#   IH02 (healthy, n=1) vs IC17/IC37/IC20/IC35 (cancer, n=4)
#
# Wilcoxon rank-sum appropriateness note:
#   With n_healthy = 1, there are exactly 5 possible rank orderings of the
#   5 observations (ranks 1–5 for the single healthy observation).  Under H₀
#   each is equally probable (p = 1/5), so the minimum achievable TWO-SIDED
#   exact p-value is 2 × (1/5) = 0.40.  The test CANNOT reach p < 0.05
#   regardless of how extreme the data are.  We run it with exact = FALSE
#   (normal approximation) for completeness — the resulting p-values are
#   displayed but are labelled as underpowered.  A valid formal test would
#   require at least n ≥ 3 in the smaller group.
# ---------------------------------------------------------------------------
message("\n--- SSP subset MDS across fragment classes ---")
message("  NOTE: Wilcoxon with n_healthy=1 — min achievable p(two-sided) = 0.40")

ssp_mds_sub <- mds_tbl |>
  filter(sample_id %in% ssp_sub_ids,
         class     %in% ssp_all_classes) |>
  mutate(
    class = factor(class, levels = ssp_all_classes),
    group = factor(if_else(sample_id == ssp_sub_ctrl, "healthy", "cancer"),
                   levels = c("healthy", "cancer"))
  )

ssp_mds_cancer  <- ssp_mds_sub |> filter(group == "cancer")
ssp_mds_healthy <- ssp_mds_sub |> filter(group == "healthy")

# Wilcoxon per class — suppress the "cannot compute exact p-value" warning
# (expected with n=1 in one group)
wilcox_ssp <- ssp_mds_sub |>
  group_by(class) |>
  summarise(
    n_cancer  = sum(group == "cancer"),
    n_healthy = sum(group == "healthy"),
    p_value   = tryCatch(
      suppressWarnings(
        wilcox.test(mds[group == "cancer"],
                    mds[group == "healthy"],
                    exact = FALSE)$p.value
      ),
      error = function(e) NA_real_
    ),
    y_top = max(mds, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    # mark everything NS explicitly — at n_healthy=1 nothing can be significant
    sig_label = if_else(
      is.na(p_value),
      "n/a",
      sprintf("p=%.2f\n(ns†)", p_value)
    )
  )

message("  Wilcoxon rank-sum (exact=FALSE, normal approx), n_cancer=4, n_healthy=1:")
for (i in seq_len(nrow(wilcox_ssp))) {
  message(sprintf("  class %-5s  p = %.4f",
                  as.character(wilcox_ssp$class[i]),
                  wilcox_ssp$p_value[i]))
}

class_xlabs <- c(
  all = "all\n(unfiltered)",
  I   = "I\n(<120 bp)",
  II  = "II\n(≥120 bp)",
  III = "III\n(35–80 bp)",
  IV  = "IV\n(120–180 bp)"
)

# y-axis: pad top to accommodate p-value text
y_range   <- range(ssp_mds_sub$mds, na.rm = TRUE)
y_pad_top <- diff(y_range) * 0.28    # ~28 % of range above highest point

p_ssp_mds <- ggplot() +
  # ---- connecting lines ----
  geom_line(data      = ssp_mds_cancer,
            aes(x     = class, y = mds,
                color = sample_id, group = sample_id),
            linewidth = 0.55, alpha = 0.65) +
  geom_line(data      = ssp_mds_healthy,
            aes(x     = class, y = mds, group = sample_id),
            linewidth = 0.8, color = group_colors["healthy"],
            alpha     = 0.8, linetype = "dashed") +
  # ---- points ----
  geom_point(data     = ssp_mds_cancer,
             aes(x    = class, y = mds, color = sample_id),
             shape    = 16, size = 3.2, alpha = 0.9) +
  geom_point(data     = ssp_mds_healthy,
             aes(x    = class, y = mds),
             shape    = 23, size = 5,
             fill     = group_colors["healthy"],
             color    = "white", stroke = 0.8) +
  # ---- Wilcoxon p-value annotation ----
  geom_text(data        = wilcox_ssp,
            aes(x       = class,
                y       = y_top + diff(y_range) * 0.04,
                label   = sig_label),
            size        = 2.6, hjust = 0.5, vjust = 0,
            color       = "grey40", lineheight = 0.9) +
  # ---- scales / labels ----
  scale_x_discrete(labels = class_xlabs) +
  scale_color_manual(values = sample_colors_sub,
                     labels = sample_labels_sub,
                     name   = "Cancer sample") +
  scale_y_continuous(
    expand = expansion(add = c(diff(y_range) * 0.04, y_pad_top))
  ) +
  labs(
    x       = "Fragment length class",
    y       = "Motif Diversity Score (MDS)",
    caption = paste0(
      "◆ ", ssp_sub_ctrl, " (healthy reference, dashed line)  |  ",
      "† Wilcoxon rank-sum (exact=FALSE); n_healthy = 1 ⇒ ",
      "min achievable p₂₋ₜₑₜ = 0.40 (underpowered)"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "right",
    plot.caption       = element_text(size = 7.5, color = "grey45", hjust = 0,
                                      margin = margin(t = 6))
  )

ggsave("analysis/ssp_subset/mds_ssp_subset.png", p_ssp_mds,
       width = 8, height = 5, dpi = 300)
message("Wrote analysis/ssp_subset/mds_ssp_subset.png")

# ---------------------------------------------------------------------------
# Per-sample log2FC scatter: class IV vs III (panel A) + class II vs I (panel B)
#   Each cancer sample vs its stratum control (same logic as lfc_dsp/ssp plots)
#     DSP           → ctrl IH01
#     SSP ≥10×      → ctrl IH02
#     SSP <10×      → ctrl IH03
#   Output: analysis/per_spl_lfc_scatters/lfc_cls_scatter_<samplename>.png
# ---------------------------------------------------------------------------
message("\n--- Per-sample log2FC class scatter plots ---")

dir.create("analysis/per_spl_lfc_scatters/IIIvsIV_IvsII", showWarnings = FALSE, recursive = TRUE)

top_label_n_scatter <- 20L   # label top N most extreme points per panel

cancer_samples <- samples |> filter(sample_group == "cancer")

walk(cancer_samples$sample_id, function(s) {
  srow <- samples |> filter(sample_id == s)
  lt   <- as.character(srow$library_type)
  cv   <- srow$coverage
  tf   <- srow$tumor_fraction
  dl   <- srow$disease_label

  stratum  <- if (lt == "DSP") "DSP" else if (cv >= 10) "SSP_ge10" else "SSP_lt10"
  ctrl_id  <- stratum_ctrl[stratum]
  lib_name <- stratum_label[stratum]

  ctrl_f <- motifs |>
    filter(sample_id == ctrl_id, as.character(class) %in% c("I","II","III","IV")) |>
    select(class, kmer, ctrl_freq = freq)

  spl_wide <- motifs |>
    filter(sample_id == s, as.character(class) %in% c("I","II","III","IV")) |>
    left_join(ctrl_f, by = c("class", "kmer")) |>
    mutate(log2fc = log2(freq / ctrl_freq)) |>
    filter(is.finite(log2fc)) |>
    select(class, kmer, first_nuc, log2fc) |>
    pivot_wider(names_from = class, values_from = log2fc, names_prefix = "lfc_") |>
    filter(!is.na(lfc_I), !is.na(lfc_II), !is.na(lfc_III), !is.na(lfc_IV))

  all_lfc <- c(spl_wide$lfc_I, spl_wide$lfc_II,
               spl_wide$lfc_III, spl_wide$lfc_IV)
  lim <- ceiling(max(abs(all_lfc), na.rm = TRUE) * 1.08 * 10) / 10

  ctrl_row      <- samples |> filter(sample_id == ctrl_id)
  ctrl_lib_name <- stratum_label[stratum]
  title_str <- sprintf("%s (%s)  %.0f×, TF≈%.0f%%  %s  vs  %s (%s)  %.0f×, TF≈%.0f%%  %s",
                       s, dl, cv, tf * 100, lib_name,
                       ctrl_id, ctrl_row$disease_label,
                       ctrl_row$coverage, ctrl_row$tumor_fraction * 100,
                       ctrl_lib_name)

  make_cls_scatter <- function(xcol, ycol, xlab, ylab) {
    d <- spl_wide |>
      mutate(diag_sum = .data[[xcol]] + .data[[ycol]])

    # Top 20 most up-shifted + top 20 most down-shifted along the diagonal
    d_label <- bind_rows(
      d |> slice_max(diag_sum, n = top_label_n_scatter, with_ties = FALSE),
      d |> slice_min(diag_sum, n = top_label_n_scatter, with_ties = FALSE)
    ) |> distinct(kmer, .keep_all = TRUE)

    d_nolabel <- d |> anti_join(d_label, by = "kmer")

    rho   <- cor(d[[xcol]], d[[ycol]], method = "spearman")
    r_lbl <- sprintf("ρ = %.2f", rho)

    ggplot(mapping = aes(x = .data[[xcol]], y = .data[[ycol]],
                         color = first_nuc)) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  linewidth = 0.3, color = "grey60") +
      geom_point(data = d_nolabel, size = 0.8, alpha = 0.35) +
      geom_point(data = d_label,   size = 2,   alpha = 0.9) +
      geom_text_repel(
        data               = d_label,
        aes(label          = as.character(kmer)),
        size               = 2.6, family = "mono",
        min.segment.length = 0, segment.size = 0.2, segment.color = "grey60",
        max.overlaps       = Inf, box.padding = 0.25
      ) +
      annotate("text", x = -lim, y = lim, label = r_lbl,
               hjust = 0, vjust = 1, size = 3.5, color = "grey25") +
      scale_color_manual(values = nuc_colors, name = "First nt") +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
      scale_x_continuous(limits = c(-lim, lim)) +
      scale_y_continuous(limits = c(-lim, lim)) +
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  }

  p_A <- make_cls_scatter(
    xcol = "lfc_IV", ycol = "lfc_III",
    xlab = "log₂ FC — class IV (120–180 bp)",
    ylab = "log₂ FC — class III (35–80 bp)"
  )
  p_B <- make_cls_scatter(
    xcol = "lfc_II", ycol = "lfc_I",
    xlab = "log₂ FC — class II (≥120 bp)",
    ylab = "log₂ FC — class I (<120 bp)"
  )

  p_out <- (p_A | p_B) +
    plot_layout(guides = "collect") +
    plot_annotation(title = title_str) &
    theme(plot.title    = element_text(size = 11, face = "plain"),
          legend.position = "right")

  outfile <- file.path("analysis/per_spl_lfc_scatters/IIIvsIV_IvsII",
                       sprintf("lfc_scatter_%s_%s_%s.png", s, dl, lt))
  ggsave(outfile, p_out, width = 13, height = 6, dpi = 300)
  message("Wrote ", outfile)
})

# ---------------------------------------------------------------------------
# Per-sample log2FC scatter: class I (x) vs class III (y) — single panel
#   Same strata/control logic and visual style as IIIvsIV_IvsII section
#   Output: analysis/per_spl_lfc_scatters/IvsIII/lfc_cls_scatter_<samplename>.png
# ---------------------------------------------------------------------------
message("\n--- Per-sample log2FC scatter: class I vs class III ---")

dir.create("analysis/per_spl_lfc_scatters/IvsIII", showWarnings = FALSE, recursive = TRUE)

walk(cancer_samples$sample_id, function(s) {
  srow <- samples |> filter(sample_id == s)
  lt   <- as.character(srow$library_type)
  cv   <- srow$coverage
  tf   <- srow$tumor_fraction
  dl   <- srow$disease_label

  stratum  <- if (lt == "DSP") "DSP" else if (cv >= 10) "SSP_ge10" else "SSP_lt10"
  ctrl_id  <- stratum_ctrl[stratum]
  lib_name <- stratum_label[stratum]

  ctrl_f <- motifs |>
    filter(sample_id == ctrl_id, as.character(class) %in% c("I", "III")) |>
    select(class, kmer, ctrl_freq = freq)

  spl_wide <- motifs |>
    filter(sample_id == s, as.character(class) %in% c("I", "III")) |>
    left_join(ctrl_f, by = c("class", "kmer")) |>
    mutate(log2fc = log2(freq / ctrl_freq)) |>
    filter(is.finite(log2fc)) |>
    select(class, kmer, first_nuc, log2fc) |>
    pivot_wider(names_from = class, values_from = log2fc, names_prefix = "lfc_") |>
    filter(!is.na(lfc_I), !is.na(lfc_III))

  lim <- ceiling(max(abs(c(spl_wide$lfc_I, spl_wide$lfc_III)), na.rm = TRUE) * 1.08 * 10) / 10

  ctrl_row      <- samples |> filter(sample_id == ctrl_id)
  title_str <- sprintf("%s (%s)  %.0f×, TF≈%.0f%%  %s  vs  %s (%s)  %.0f×, TF≈%.0f%%  %s",
                       s, dl, cv, tf * 100, lib_name,
                       ctrl_id, ctrl_row$disease_label,
                       ctrl_row$coverage, ctrl_row$tumor_fraction * 100,
                       lib_name)

  d <- spl_wide |> mutate(diag_sum = lfc_I + lfc_III)
  d_label <- bind_rows(
    d |> slice_max(diag_sum, n = top_label_n_scatter, with_ties = FALSE),
    d |> slice_min(diag_sum, n = top_label_n_scatter, with_ties = FALSE)
  ) |> distinct(kmer, .keep_all = TRUE)
  d_nolabel <- d |> anti_join(d_label, by = "kmer")

  rho   <- cor(d$lfc_I, d$lfc_III, method = "spearman")
  r_lbl <- sprintf("ρ = %.2f", rho)

  p_out <- ggplot(mapping = aes(x = lfc_I, y = lfc_III, color = first_nuc)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, color = "grey60") +
    geom_point(data = d_nolabel, size = 0.8, alpha = 0.35) +
    geom_point(data = d_label,   size = 2,   alpha = 0.9) +
    geom_text_repel(
      data               = d_label,
      aes(label          = as.character(kmer)),
      size               = 2.6, family = "mono",
      min.segment.length = 0, segment.size = 0.2, segment.color = "grey60",
      max.overlaps       = Inf, box.padding = 0.25
    ) +
    annotate("text", x = -lim, y = lim, label = r_lbl,
             hjust = 0, vjust = 1, size = 3.5, color = "grey25") +
    scale_color_manual(values = nuc_colors, name = "First nt") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    scale_x_continuous(limits = c(-lim, lim)) +
    scale_y_continuous(limits = c(-lim, lim)) +
    labs(title = title_str,
         x     = "log₂ FC — class I (<120 bp)",
         y     = "log₂ FC — class III (35–80 bp)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor  = element_blank(),
          plot.title        = element_text(size = 11, face = "plain"),
          legend.position   = "right")

  outfile <- file.path("analysis/per_spl_lfc_scatters/IvsIII",
                       sprintf("lfc_scatter_%s_%s_%s.png", s, dl, lt))
  ggsave(outfile, p_out, width = 7, height = 6, dpi = 300)
  message("Wrote ", outfile)
})
