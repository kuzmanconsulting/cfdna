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

ggsave("analysis/coverage_vs_tf.png", p_cov_tf,
       width = 5, height = 4, dpi = 300)
message("\nWrote analysis/coverage_vs_tf.png")

# ---------------------------------------------------------------------------
# Fragment length ridgeplot (35-500 bp)
# ---------------------------------------------------------------------------
cache_dir  <- "analysis/.cache"
cache_file <- file.path(cache_dir, "frag_lengths_35_500.tsv")

if (!file.exists(cache_file)) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
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

ggsave("analysis/frag_length_ridges.png", p_frag,
       width = 8, height = 0.35 * nrow(samples) + 2, dpi = 300)
message("Wrote analysis/frag_length_ridges.png")

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

ggsave("analysis/motif_freq_ge10.png", p_top_mot,
       width = 12, height = 7, dpi = 300)
message("Wrote analysis/motif_freq_ge10.png")

# Bottom: <10x, faceted by group only (SSP only — no DSP samples below 10×)
bot_d <- motifs |> filter(coverage < 10)
stopifnot(all(bot_d$library_type == "SSP"))

p_bot_mot <- make_motif_bar(bot_d, "sample_group") +
  ggtitle("Median 4-mer end-motif frequencies, coverage < 10× (SSP only)")

ggsave("analysis/motif_freq_lt10.png", p_bot_mot,
       width = 7, height = 7, dpi = 300)
message("Wrote analysis/motif_freq_lt10.png")

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

ggsave("analysis/mds_strip.png", p_mds,
       width = 8, height = 4, dpi = 300)
message("Wrote analysis/mds_strip.png")

# Per-sample plots
if (isTRUE(config$skip_per_spl_4mer_freq)) {
  message("Skipping per-sample 4-mer plots (config$skip_per_spl_4mer_freq = TRUE)")
} else {
per_sample_dir <- "analysis/4mer_freqs"
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
  ggsave(file.path("analysis", out_files[st]), p, width = 8, height = 7, dpi = 300)
  message("Wrote analysis/", out_files[st])
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

ggsave("analysis/pca_all_samples.png", p_pca_all, width = 7, height = 6, dpi = 300)
message("Wrote analysis/pca_all_samples.png")

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

ggsave("analysis/pca_by_library.png", p_pca_libs, width = 14, height = 6, dpi = 300)
message("Wrote analysis/pca_by_library.png")
