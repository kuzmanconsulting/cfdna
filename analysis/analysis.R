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
                              disease_label),
            by = "sample_id") |>
  mutate(label = if_else(sample_group == "healthy",
                         sample_id, paste(sample_id, disease_label)))

mds_class <- mds_tbl |>
  filter(class %in% classes) |>
  mutate(class = factor(class, levels = classes),
         group_key = interaction(sample_group, library_type, drop = TRUE))

dodge_w <- 0.7
g_lvls  <- levels(mds_class$group_key)
n_g     <- length(g_lvls)
offsets <- setNames(-dodge_w/2 + (seq_len(n_g) - 0.5) * dodge_w / n_g, g_lvls)

mds_class <- mds_class |>
  mutate(x_num = as.integer(class) + offsets[as.character(group_key)])

mds_all_ref <- mds_tbl |> filter(class == "all")

p_mds <- ggplot(mds_class, aes(x_num, mds,
                               color = sample_group, shape = library_type)) +
  geom_hline(data = mds_all_ref,
             aes(yintercept = mds, color = sample_group),
             linetype = "dashed", linewidth = 0.3, alpha = 0.4,
             show.legend = FALSE) +
  geom_point(size = 2, stroke = 0.8) +
  geom_text_repel(aes(label = label),
                  size = 2, color = "grey40", alpha = 0.8,
                  min.segment.length = 0, segment.size = 0.2,
                  segment.color = "grey70", max.overlaps = Inf,
                  box.padding = 0.3) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = c(SSP = 16, DSP = 17)) +
  scale_x_continuous(breaks = seq_along(classes), labels = classes) +
  labs(x = "Length class", y = "Motif diversity score",
       color = "Group", shape = "Library") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("analysis/mds_strip.png", p_mds,
       width = 7, height = 5, dpi = 300)
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
