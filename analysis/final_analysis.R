#!/usr/bin/env Rscript
# Final end-motif analysis: 5 samples (1 healthy + 4 cancer), 4 length classes.
# Outputs are numbered by analytical question:
#   Q1 — do 4-mers differ between cancer and normal? (classes I-IV)
#     Q1_1_MDS_strip.png — per-sample MDS by class (sanity: all near ~0.95)
#
# Inputs  (read-only): cfdna-finale-snakemake/results/<S>/end_motifs/
#   <S>.<cls>.tsv      headerless: kmer <tab> freq  (256 rows, sums to 1)
#   <S>.<cls>.mds.txt  single float (normalized Shannon entropy)
# Outputs: analysis/Q*.{png,tsv}

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidytext)
})

# --- paths -------------------------------------------------------------------
here    <- normalizePath(".")
res_dir <- file.path(here, "cfdna-finale-snakemake", "results")
if (!dir.exists(res_dir))
  res_dir <- file.path(here, "..", "cfdna-finale-snakemake", "results")
stopifnot(dir.exists(res_dir))
out_dir <- if (dir.exists(file.path(here, "analysis"))) file.path(here, "analysis") else here

# --- sample registry ---------------------------------------------------------
ss_path <- file.path(out_dir, "samplesheet.csv")
if (!file.exists(ss_path)) stop("samplesheet.csv not found in analysis/: ", ss_path)
samples <- read_csv(ss_path, col_types = "cccc")
classes <- c("I", "II", "III", "IV")

sample_colors <- c(
  IH02 = "#2166AC",
  IC37 = "#B2182B", IC33 = "#D6604D",
  IC15 = "#1B7837", IC32 = "#7FBC41"
)

# --- load MDS ----------------------------------------------------------------
read_mds <- function(s, cls) {
  f <- file.path(res_dir, s, "end_motifs", sprintf("%s.%s.mds.txt", s, cls))
  if (!file.exists(f)) { warning("missing: ", f); return(NA_real_) }
  as.numeric(readLines(f, n = 1))
}

mds <- tidyr::expand_grid(sample = samples$sample, class = classes) |>
  mutate(mds = map2_dbl(sample, class, read_mds)) |>
  left_join(samples, by = "sample")

mds_all <- tibble(sample = samples$sample,
                  mds    = map_dbl(samples$sample, \(s) read_mds(s, "all"))) |>
  left_join(samples, by = "sample")

# --- load 4-mer frequencies --------------------------------------------------
read_motifs <- function(s, cls) {
  f <- file.path(res_dir, s, "end_motifs", sprintf("%s.%s.tsv", s, cls))
  if (!file.exists(f)) { warning("missing: ", f); return(NULL) }
  read_tsv(f, col_names = c("kmer", "freq"), col_types = "cd") |>
    mutate(sample = s, class = cls)
}

motifs <- tidyr::expand_grid(sample = samples$sample, class = classes) |>
  pmap(\(sample, class) read_motifs(sample, class)) |>
  list_rbind() |>
  left_join(samples, by = "sample")

# ============================================================================
# Q1_1 — per-sample MDS by class
# ============================================================================
p_mds <- mds |>
  mutate(class = factor(class, levels = classes),
         sample = factor(sample, levels = samples$sample)) |>
  ggplot(aes(class, mds, color = sample, shape = group)) +
  geom_hline(data = mds_all, aes(yintercept = mds, color = sample),
             linetype = "dashed", linewidth = 0.3, show.legend = FALSE) +
  geom_point(size = 2, stroke = 0.8,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = sample_colors) +
  scale_shape_manual(values = c(healthy = 1, cancer = 2)) +
  scale_y_continuous(limits = c(0.8, 1), breaks = c(0, 0.8, 0.9, 0.925, 0.95, 0.975, 1)) +
  labs(title    = "Normalized Shannon entropy by length class",
       subtitle = "(dashed lines - all fragments)",
       x = "length class", y = "Motif diversity score",
       color = "Sample", shape = "Group") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Q1_1_MDS_strip.png"), p_mds,
       width = 5, height = 4, dpi = 300)

message("Wrote Q1_1_MDS_strip.png to ", normalizePath(out_dir))

# ============================================================================
# Q1_2 — Jensen-Shannon divergence: each cancer vs healthy IH02, per class
# ============================================================================
js_div <- function(p, q) {
  p <- p / sum(p); q <- q / sum(q); m <- 0.5 * (p + q)
  kl <- function(a, b) sum(ifelse(a > 0, a * log2(a / b), 0))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

vec <- \(s, cls) motifs |>
  filter(sample == s, class == cls) |>
  arrange(kmer) |>
  pull(freq)

cancers <- samples$sample[samples$group == "cancer"]

read_n_ends <- function(s, cls) {
  f <- file.path(res_dir, s, "end_motifs", sprintf("%s.%s.n_ends.txt", s, cls))
  if (!file.exists(f)) { warning("missing: ", f); return(NA_real_) }
  as.numeric(readLines(f, n = 1))
}

# Analytic JS noise floor (bits): expected JS between two independent multinomial
# draws (sizes na, nb over K=256 bins) from a COMMON distribution.
js_floor <- \(na, nb, K = 256) (K - 1) / (4 * log(2)) * (1 / na + 1 / nb)

jsd <- tidyr::expand_grid(class = classes, sample = cancers) |>
  mutate(
    JSD            = pmap_dbl(list(sample, class), \(s, cls) js_div(vec(s, cls), vec("IH02", cls))),
    n_cancer       = map2_dbl(sample, class, read_n_ends),
    n_healthy      = map_dbl(class, \(cls) read_n_ends("IH02", cls)),
    JS_noise_floor = js_floor(n_cancer, n_healthy),
    JS_over_floor  = JSD / JS_noise_floor,
  ) |>
  left_join(samples, by = "sample")

write_tsv(jsd, file.path(out_dir, "Q1_2_JSD.tsv"))

p_jsd <- jsd |>
  mutate(class  = factor(class, levels = classes),
         sample = factor(sample, levels = cancers)) |>
  ggplot(aes(class, JSD, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = sample_colors[cancers]) +
  labs(title    = "End-motif divergence from healthy (IH02)",
       subtitle = "Jensen-Shannon divergence, per length class",
       x = "length class", y = "JS divergence (bits)",
       fill = "Cancer sample") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Q1_2_JSD_bar.png"), p_jsd,
       width = 6, height = 4, dpi = 300)

message("Wrote Q1_2_JSD_bar.png to ", normalizePath(out_dir))

# ============================================================================
# Q1_3 — JSD normalized by analytic noise floor (signal-to-noise ratio)
# ============================================================================
p_jsd_norm <- jsd |>
  mutate(class  = factor(class, levels = classes),
         sample = factor(sample, levels = cancers)) |>
  ggplot(aes(class, JS_over_floor, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = sample_colors[cancers]) +
  scale_y_continuous(trans = "log10") +
  labs(title    = "JSD normalized by analytic sampling-noise floor",
       x = "length class", y = "JSD / noise floor (log10)",
       fill = "Cancer sample") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Q1_3_JSD_normalized_bar.png"), p_jsd_norm,
       width = 6, height = 4, dpi = 300)

message("Wrote Q1_3_JSD_normalized_bar.png to ", normalizePath(out_dir))

# ============================================================================
# Q1_4 — PCA biplot per length class (samples = points, top loadings = arrows)
# ============================================================================
top_n_arr <- 15L

pca_data <- map(classes, function(cls) {
  wide <- motifs |>
    filter(class == cls) |>
    select(sample, kmer, freq) |>
    pivot_wider(names_from = kmer, values_from = freq)

  mat <- wide |> select(-sample) |> as.matrix()
  rownames(mat) <- wide$sample

  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  scores <- as_tibble(pca$x[, 1:2], rownames = "sample") |>
    left_join(samples, by = "sample") |>
    mutate(class = cls)

  load_scale <- max(abs(pca$x[, 1:2])) / max(abs(pca$rotation[, 1:2]))
  top_kmers  <- order(sqrt(pca$rotation[, 1]^2 + pca$rotation[, 2]^2),
                      decreasing = TRUE)[seq_len(top_n_arr)]
  loadings <- as_tibble(pca$rotation[top_kmers, 1:2], rownames = "kmer") |>
    mutate(across(c(PC1, PC2), \(x) x * load_scale), class = cls)

  list(scores = scores, loadings = loadings,
       pc1_var = pct[1], pc2_var = pct[2])
})

scores_all   <- map(pca_data, "scores")   |> list_rbind()
loadings_all <- map(pca_data, "loadings") |> list_rbind()

p_pca <- ggplot() +
  geom_segment(data = loadings_all,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.12, "cm")),
               linewidth = 0.3, color = "grey60", alpha = 0.7) +
  geom_text(data = loadings_all,
            aes(PC1, PC2, label = kmer),
            size = 2, color = "grey40", vjust = -0.4) +
  geom_point(data = scores_all,
             aes(PC1, PC2, color = sample, shape = group), size = 3) +
  geom_text(data = scores_all,
            aes(PC1, PC2, label = sample, color = sample),
            nudge_y = 0.003, size = 3, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = sample_colors) +
  scale_shape_manual(values = c(healthy = 16, cancer = 17)) +
  facet_wrap(~ factor(class, levels = classes), scales = "free", ncol = 2) +
  labs(title    = "End-motif PCA biplot by length class",
       subtitle = sprintf("Top %d 4-mer loadings shown per panel  |  n = 5 samples", top_n_arr),
       x = "PC1", y = "PC2", color = "Sample", shape = "Group") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Q1_4_PCA_biplot.png"), p_pca,
       width = 10, height = 8, dpi = 300)

message("Wrote Q1_4_PCA_biplot.png to ", normalizePath(out_dir))

# ============================================================================
# Q1_5 — per-motif log2 fold-change: each cancer vs healthy IH02
# Top N motifs selected per class by max |log2FC| across the 4 cancers.
# ============================================================================
top_n_mot <- 20L
eps       <- 1e-9

log2fc <- motifs |>
  select(sample, class, kmer, freq) |>
  pivot_wider(names_from = sample, values_from = freq) |>
  mutate(across(all_of(cancers), \(x) log2((x + eps) / (IH02 + eps)),
                .names = "lfc_{.col}")) |>
  select(class, kmer, starts_with("lfc_")) |>
  pivot_longer(starts_with("lfc_"),
               names_to = "sample", names_prefix = "lfc_",
               values_to = "log2FC")

write_tsv(log2fc, file.path(out_dir, "Q1_5_log2FC.tsv"))

top_motifs <- log2fc |>
  group_by(class, kmer) |>
  summarise(max_abs = max(abs(log2FC)), mean_lfc = mean(log2FC), .groups = "drop") |>
  group_by(class) |>
  slice_max(max_abs, n = top_n_mot) |>
  select(class, kmer, mean_lfc)

plot_df <- log2fc |>
  inner_join(top_motifs, by = c("class", "kmer")) |>
  mutate(class  = factor(class, levels = classes),
         sample = factor(sample, levels = cancers))

p_lfc <- ggplot(plot_df, aes(log2FC, reorder_within(kmer, mean_lfc, class),
                             fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = sample_colors[cancers]) +
  tidytext::scale_y_reordered() +
  facet_wrap(~ class, scales = "free_y", ncol = 2) +
  labs(title = sprintf("Top %d 4-mers by |log2FC| vs healthy, per length class", top_n_mot),
       x = "log2(cancer / IH02)", y = "4-mer",
       fill = "Cancer sample") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(family = "mono", size = 7))

ggsave(file.path(out_dir, "Q1_5_log2FC_top.png"), p_lfc,
       width = 6, height = 7, dpi = 300)

message("Wrote Q1_5_log2FC_top.png to ", normalizePath(out_dir))

# ============================================================================
# Q1_6 — Pan-cancer motif consistency: direction agreement across 4 cancers
# Stacked bar: per class, how many motifs have 4/4, 3/4, 2/4, 1/4, 0/4
# cancer samples agreeing on enrichment (positive log2FC) vs depletion.
# ============================================================================
consistency <- log2fc |>
  group_by(class, kmer) |>
  summarise(
    n_up   = sum(log2FC > 0),
    n_down = sum(log2FC < 0),
    .groups = "drop"
  ) |>
  mutate(
    tier = case_when(
      n_up  == 4 ~ "4/4 enriched",
      n_down == 4 ~ "4/4 depleted",
      n_up  == 3 ~ "3/4 enriched",
      n_down == 3 ~ "3/4 depleted",
      TRUE        ~ "mixed (≤2/4)"
    ),
    tier = factor(tier, levels = c("4/4 enriched", "3/4 enriched",
                                   "mixed (≤2/4)",
                                   "3/4 depleted", "4/4 depleted"))
  )

tier_colors <- c(
  "4/4 enriched" = "#B2182B",
  "3/4 enriched" = "#D6604D",
  "mixed (≤2/4)" = "grey80",
  "3/4 depleted" = "#4393C3",
  "4/4 depleted" = "#2166AC"
)

p_cons <- consistency |>
  count(class, tier) |>
  mutate(class = factor(class, levels = classes)) |>
  ggplot(aes(class, n, fill = tier)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tier_colors) +
  labs(title = "Pan-cancer 4-mer consistency by length class",
       x = "length class", y = "number of 4-mers (of 256)",
       fill = "Consistency") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Q1_6_consistency_bar.png"), p_cons,
       width = 4, height = 4, dpi = 300)

write_tsv(consistency, file.path(out_dir, "Q1_6_consistency.tsv"))

message("Wrote Q1_6_consistency_bar.png to ", normalizePath(out_dir))

# ============================================================================
# Q1_6a — motifs with |mean log2FC| > threshold, colored by consistency tier,
#          RC partner annotated on x-axis label.
# ============================================================================
rc <- function(kmer) {
  comp <- chartr("ACGT", "TGCA", kmer)
  paste(rev(strsplit(comp, "")[[1]]), collapse = "")
}

kmer_mean_lfc <- log2fc |>
  group_by(class, kmer) |>
  summarise(mean_lfc = mean(log2FC), .groups = "drop")

lfc_thresh <- 0.4

cons_lfc <- consistency |>
  left_join(kmer_mean_lfc, by = c("class", "kmer")) |>
  filter(abs(mean_lfc) > lfc_thresh) |>
  mutate(class = factor(class, levels = classes)) |>
  mutate(
    rc_kmer  = map_chr(kmer, rc),
    rc_status = pmap_chr(list(rc_kmer, class, mean_lfc), \(rk, cls, lfc) {
      partner <- kmer_mean_lfc |> filter(kmer == rk, class == cls)
      if (nrow(partner) == 0 || abs(partner$mean_lfc) <= lfc_thresh)
        return("absent")
      if (sign(partner$mean_lfc) == sign(lfc)) "same" else "opposite"
    }),
    # x-axis label: kmer with RC partner appended if present in plot
    x_label = if_else(rc_status != "absent",
                      paste0(kmer, "\n↔", rc_kmer),
                      kmer)
  )

rc_border <- c(same = "#01ff3c", opposite = "#FF00CC", absent = NA)

p_cons_bar <- cons_lfc |>
  ggplot(aes(reorder_within(kmer, mean_lfc, class), mean_lfc, fill = tier)) +
  geom_col(aes(color = rc_status), width = 0.8, linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_fill_manual(values = tier_colors) +
  scale_color_manual(values = rc_border, na.value = NA,
                     breaks = c("same", "opposite"),
                     labels = c(same = "RC also shifted (same dir)",
                                opposite = "RC also shifted (opposite)"),
                     name = "RC partner",
                     guide = guide_legend(override.aes = list(fill = "white", linewidth = 1))) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  scale_x_reordered() +
  facet_wrap(~ class, scales = "free_x", ncol = 2) +
  labs(title = sprintf("4-mers with |mean log2FC| > %.1f vs healthy", lfc_thresh),
       x = "4-mer",
       y = "mean log2FC (cancer / IH02)",
       fill = "Consistency") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(family = "mono", size = 7, angle = 45, hjust = 1))

ggsave(file.path(out_dir, "Q1_6a_consistency_motifs.png"), p_cons_bar,
       width = 12, height = 5.5, dpi = 300)

message("Wrote Q1_6a_consistency_motifs.png to ", normalizePath(out_dir))
