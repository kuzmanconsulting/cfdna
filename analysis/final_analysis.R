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
  library(ggrepel)
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
# Q1_0 — per-sample 4-mer frequency barplot (fill = first nucleotide)
# ============================================================================
nuc_colors <- c(A = "#4393C3", T = "#FFC107", C = "#D6604D", G = "#4DAC26")

motif_bar_data <- motifs |>
  mutate(
    first_nuc = substr(kmer, 1, 1),
    kmer      = factor(kmer, levels = sort(unique(kmer)))
  )

p_motif_bar <- motif_bar_data |>
  mutate(sample = factor(sample, levels = samples$sample)) |>
  ggplot(aes(x = kmer, y = freq, fill = first_nuc)) +
  geom_col(width = 1, linewidth = 0) +
  scale_fill_manual(values = nuc_colors, name = "First nt") +
  scale_x_discrete(breaks = NULL) +
  facet_grid(class ~ sample) +
  labs(
    title = "4-mer end-motif frequencies",
    x     = "4-mer motif (256)",
    y     = "Frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )

ggsave(file.path(out_dir, "Q1_0_motif_freq_bar.png"), p_motif_bar,
       width = 3 * length(unique(motif_bar_data$sample)), height = 6, dpi = 300)
message("Wrote Q1_0_motif_freq_bar.png to ", normalizePath(out_dir))

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

# ============================================================================
# Q1_7 — PCA on cross-class consistently shifted motifs
# Feature set: motifs that are 4/4 enriched OR 4/4 depleted in ALL 4 classes.
# Observations: all 5 samples × 4 classes (20 rows).
# Samples are labeled; classes are indicated by point shape / facet stripe.
# ============================================================================

# Per-class consistent motifs: 4/4 enriched or 4/4 depleted within that class
pca7_data <- map(classes, function(cls) {
  cons_cls <- consistency |>
    filter(class == cls, tier %in% c("4/4 enriched", "4/4 depleted")) |>
    select(kmer, direction = tier) |>
    mutate(direction = if_else(direction == "4/4 enriched", "enriched", "depleted"))

  message(sprintf("Q1_7 class %s: %d consistent motifs (%d enriched, %d depleted)",
                  cls, nrow(cons_cls),
                  sum(cons_cls$direction == "enriched"),
                  sum(cons_cls$direction == "depleted")))

  if (nrow(cons_cls) < 2) {
    message("  → skipping (fewer than 2 motifs)")
    return(NULL)
  }

  feat_wide <- motifs |>
    filter(class == cls, kmer %in% cons_cls$kmer) |>
    select(sample, kmer, freq) |>
    pivot_wider(names_from = kmer, values_from = freq)

  mat <- feat_wide |> select(-sample) |> as.matrix()
  rownames(mat) <- feat_wide$sample

  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  scores <- as_tibble(pca$x[, 1:2], rownames = "sample") |>
    left_join(samples, by = "sample") |>
    mutate(class = cls)

  load_scale <- max(abs(pca$x[, 1:2])) / max(abs(pca$rotation[, 1:2]))
  top_idx <- order(sqrt(pca$rotation[, 1]^2 + pca$rotation[, 2]^2),
                   decreasing = TRUE)[seq_len(min(20L, nrow(pca$rotation)))]
  loadings <- as_tibble(pca$rotation[top_idx, 1:2], rownames = "kmer") |>
    mutate(across(c(PC1, PC2), \(x) x * load_scale)) |>
    left_join(cons_cls, by = "kmer") |>
    mutate(class = cls)

  list(scores = scores, loadings = loadings,
       pc1_var = pct[1], pc2_var = pct[2], cls = cls)
})
names(pca7_data) <- classes
pca7_data <- Filter(Negate(is.null), pca7_data)

scores7_all   <- map(pca7_data, "scores")   |> list_rbind()
loadings7_all <- map(pca7_data, "loadings") |> list_rbind()

var_labels <- map_dfr(pca7_data, \(d) tibble(
  class   = d$cls,
  pc1_var = d$pc1_var,
  pc2_var = d$pc2_var
))

# per-panel axis labels via dummy aesthetics on an empty geom
axis_df <- scores7_all |>
  distinct(class) |>
  left_join(var_labels, by = "class") |>
  mutate(x_lab = sprintf("PC1 (%.1f%%)", pc1_var),
         y_lab = sprintf("PC2 (%.1f%%)", pc2_var))

nudge_y7 <- scores7_all |>
  group_by(class) |>
  summarise(ny = diff(range(PC2)) * 0.08, .groups = "drop")

scores7_all <- scores7_all |> left_join(nudge_y7, by = "class")

p_q17 <- ggplot() +
  geom_segment(data = loadings7_all,
               aes(x = 0, y = 0, xend = PC1, yend = PC2, color = direction),
               arrow = arrow(length = unit(0.12, "cm")),
               linewidth = 0.4, alpha = 0.7) +
  geom_text(data = loadings7_all,
            aes(PC1, PC2, label = kmer, color = direction),
            size = 2, vjust = -0.4, show.legend = FALSE) +
  geom_point(data = scores7_all,
             aes(PC1, PC2, color = sample, shape = group), size = 3) +
  geom_text(data = scores7_all,
            aes(PC1, PC2, label = sample, color = sample, nudge_y = ny),
            stat = "identity", vjust = -0.6,
            size = 2.8, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(
    values = c(sample_colors, enriched = "#B2182B", depleted = "#2166AC"),
    breaks = names(sample_colors), name = "Sample") +
  scale_shape_manual(values = c(healthy = 16, cancer = 17), name = "Group") +
  facet_wrap(~ factor(class, levels = classes), scales = "free", ncol = 2,
             labeller = as_labeller(\(x) {
               d <- var_labels[var_labels$class == x, ]
               sprintf("Class %s  |  PC1 %.1f%%  PC2 %.1f%%", x, d$pc1_var, d$pc2_var)
             })) +
  labs(title = "PCA on 4/4-consistent 4-mers, per length class",
       x = "PC1", y = "PC2") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Q1_7_PCA_consistent_motifs.png"), p_q17,
       width = 10, height = 8, dpi = 300)

message("Wrote Q1_7_PCA_consistent_motifs.png to ", normalizePath(out_dir))

# ============================================================================
# Q2 — Does coverage depth distort motif profiles?
#   Within each cancer subtype, compare deep vs shallow sample motif profiles.
#   If depth noise << cancer-vs-normal signal, shallow samples are trustworthy.
#
# Requires samplesheet columns: subtype, depth ("deep"/"shallow")
# One deep + one shallow representative per subtype; if multiple exist the
# script averages their frequencies (future-proof for larger sample lists).
#
# Outputs:
#   Q2_1_depth_scatter.png  — 256-kmer scatter: deep (x) vs shallow (y) per
#                             subtype × class; Pearson r annotated
#   Q2_2_depth_vs_signal.png — JSD bar: depth-noise vs cancer-signal per
#                              subtype × class
# ============================================================================

# --- helpers -----------------------------------------------------------------
jsd <- function(p, q) {
  m <- (p + q) / 2
  kl <- function(a, b) sum(ifelse(a > 0, a * log2(a / b), 0))
  (kl(p, m) + kl(q, m)) / 2
}

avg_profile <- function(samples_vec, cls) {
  freqs <- map(samples_vec, \(s) read_motifs(s, cls))
  freqs <- keep(freqs, \(x) !is.null(x))
  if (length(freqs) == 0) return(NULL)
  bind_rows(freqs) |>
    group_by(kmer) |>
    summarise(freq = mean(freq), .groups = "drop")
}

# --- build depth-pair data ---------------------------------------------------
subtypes_q2 <- samples |>
  filter(group == "cancer", !is.na(subtype), subtype != "—") |>
  distinct(subtype) |>
  pull(subtype)

depth_pairs <- tidyr::expand_grid(subtype = subtypes_q2, class = classes) |>
  mutate(
    deep_samples    = map(subtype, \(st) samples |>
                            filter(group == "cancer", subtype == st, depth == "deep") |>
                            pull(sample)),
    shallow_samples = map(subtype, \(st) samples |>
                            filter(group == "cancer", subtype == st, depth == "shallow") |>
                            pull(sample)),
    deep_profile    = map2(deep_samples,    class, avg_profile),
    shallow_profile = map2(shallow_samples, class, avg_profile)
  ) |>
  filter(!map_lgl(deep_profile, is.null), !map_lgl(shallow_profile, is.null))

scatter_data <- depth_pairs |>
  mutate(
    joined = map2(deep_profile, shallow_profile, \(d, s)
      inner_join(d, s, by = "kmer", suffix = c("_deep", "_shallow"))
    )
  ) |>
  select(subtype, class, joined) |>
  unnest(joined)

# Pearson r per panel
r_labels <- scatter_data |>
  group_by(subtype, class) |>
  summarise(r = cor(freq_deep, freq_shallow), .groups = "drop") |>
  mutate(label = sprintf("r = %.3f", r))

# axis range (shared across all panels for comparability)
ax_max <- max(c(scatter_data$freq_deep, scatter_data$freq_shallow)) * 1.05

# GC content for point colouring (number of C/G bases in 4-mer)
scatter_data <- scatter_data |>
  mutate(gc = str_count(kmer, "[GC]"))

high_motifs <- scatter_data |> filter(freq_deep > 0.015 | freq_shallow > 0.015)

p_scatter <- ggplot(scatter_data, aes(x = freq_deep, y = freq_shallow, colour = factor(gc))) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 1.4, alpha = 0.8) +
  geom_label_repel(data = high_motifs, aes(x = freq_deep, y = freq_shallow, label = kmer),
                   size = 2.5, colour = "black", fill = alpha("white", 0.8),
                   label.size = 0.2, min.segment.length = 0,
                   box.padding = 0.3, inherit.aes = FALSE) +
  geom_text(data = r_labels, aes(label = label),
            x = ax_max * 0.05, y = ax_max * 0.95,
            hjust = 0, vjust = 1, size = 3, colour = "black", inherit.aes = FALSE) +
  scale_colour_manual(
    values = c("0" = "#4575B4", "1" = "#74ADD1", "2" = "#FEE090",
               "3" = "#F46D43", "4" = "#D73027"),
    name = "GC count"
  ) +
  scale_x_continuous(limits = c(0, ax_max), labels = scales::label_number(accuracy = 0.001)) +
  scale_y_continuous(limits = c(0, ax_max), labels = scales::label_number(accuracy = 0.001)) +
  facet_grid(subtype ~ class,
             labeller = labeller(class = \(x) paste("Class", x))) +
  labs(
    title  = "Coverage depth: deep vs shallow motif profile concordance",
    x      = "Deep sample — 4-mer frequency",
    y      = "Shallow sample — 4-mer frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey92"),
    legend.position   = "right"
  )

ggsave(file.path(out_dir, "Q2_1_depth_scatter.png"), p_scatter,
       width = 10, height = 3 * length(subtypes_q2), dpi = 300)
message("Wrote Q2_1_depth_scatter.png to ", normalizePath(out_dir))

# --- Q2_2: depth-noise JSD vs cancer-signal JSD ------------------------------
# depth noise  = JSD(deep_cancer, shallow_cancer) per subtype × class
# cancer signal = JSD(deep_cancer, IH02) per subtype × class

ih02_profiles <- map(classes, \(cls) {
  p <- read_motifs("IH02", cls)
  if (is.null(p)) return(NULL)
  tibble(class = cls, kmer = p$kmer, freq_ref = p$freq)
}) |> list_rbind()

jsd_bars <- depth_pairs |>
  mutate(
    jsd_depth = map2_dbl(deep_profile, shallow_profile, \(d, s) {
      j <- inner_join(d, s, by = "kmer")
      jsd(j$freq.x, j$freq.y)
    }),
    jsd_signal = pmap_dbl(list(deep_profile, class), \(d, cls) {
      ref <- ih02_profiles |> filter(class == cls) |> select(kmer, freq_ref)
      if (nrow(ref) == 0) return(NA_real_)
      j <- inner_join(d, ref, by = "kmer")
      jsd(j$freq, j$freq_ref)
    })
  ) |>
  select(subtype, class, jsd_depth, jsd_signal) |>
  pivot_longer(c(jsd_depth, jsd_signal),
               names_to = "comparison", values_to = "jsd") |>
  mutate(
    comparison = factor(comparison,
      levels = c("jsd_depth", "jsd_signal"),
      labels = c("Depth noise\n(deep vs shallow)", "Cancer signal\n(deep vs IH02)")
    )
  )

p_jsd2 <- ggplot(jsd_bars, aes(x = class, y = jsd, fill = comparison)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(
    values = c("Depth noise\n(deep vs shallow)" = "#BDBDBD",
               "Cancer signal\n(deep vs IH02)"  = "#B2182B"),
    name = NULL
  ) +
  facet_wrap(~subtype, ncol = length(subtypes_q2)) +
  labs(
    title = "Depth noise vs cancer signal (JSD)",
    x     = "Length class",
    y     = "Jensen-Shannon divergence"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "right"
  )

ggsave(file.path(out_dir, "Q2_2_depth_vs_signal.png"), p_jsd2,
       width = 4 * length(subtypes_q2), height = 4, dpi = 300)
message("Wrote Q2_2_depth_vs_signal.png to ", normalizePath(out_dir))
