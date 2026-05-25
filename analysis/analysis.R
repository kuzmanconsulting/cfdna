#!/usr/bin/env Rscript
# Load samples from cfdna-finale-snakemake/samplesheet.csv, attach disease
# (from snyder2016_metadata_GSE.csv) and disease_label (from
# analysis/tumor_types.csv).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggridges)
})

samples      <- read_csv("cfdna-finale-snakemake/samplesheet.csv", show_col_types = FALSE)
meta         <- read_csv("snyder2016_metadata_GSE.csv",            show_col_types = FALSE)
tumor_types  <- read_csv("analysis/tumor_types.csv",               show_col_types = FALSE)

samples <- samples |>
  left_join(meta |> select(sample_id, disease), by = "sample_id") |>
  left_join(tumor_types, by = "disease") |>
  rename(disease_label = label)

missing <- samples |> filter(is.na(disease_label))
if (nrow(missing))
  stop("No disease_label for: ", paste(unique(missing$disease), collapse = "; "))

print(samples, n = Inf)

# Coverage vs tumor fraction overview
p_cov_tf <- ggplot(samples, aes(coverage, tumor_fraction,
                                shape = library_type, color = sample_group)) +
  geom_point(size = 2.2) +
  geom_text_repel(aes(label = if_else(sample_group == "healthy",
                                       sample_id,
                                       paste(sample_id, disease_label))),
                  size = 2.6, color = "grey40", alpha = 0.7,
                  min.segment.length = 0, segment.size = 0.2,
                  segment.color = "grey70", max.overlaps = Inf) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Coverage (log10)", y = "Tumor fraction (ichorCNA)",
       shape = "Library", color = "Group") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("analysis/coverage_vs_tf.png", p_cov_tf,
       width = 5, height = 4.5, dpi = 300)
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
  group_by(sample_id) |>
  mutate(height = count / max(count)) |>   # free-y: normalize each ridge to its own max
  ungroup()

sample_order <- samples |>
  mutate(group_rank = if_else(sample_group == "healthy", 0L, 1L)) |>
  arrange(group_rank, desc(coverage)) |>
  mutate(y_label = paste0(sample_id, " (", disease_label, ")")) |>
  select(sample_id, y_label)

frag <- frag |>
  left_join(sample_order, by = "sample_id") |>
  mutate(y_label = factor(y_label, levels = rev(sample_order$y_label)))

lib_colors <- c(SSP = "#2166AC", DSP = "#B2182B")

p_frag <- ggplot(frag, aes(x = length, y = y_label, height = height,
                            fill = library_type, color = library_type)) +
  annotate("rect", xmin = 35,  xmax = 80,  ymin = -Inf, ymax = Inf,
           fill = "#FDD0A2", alpha = 0.35) +
  annotate("rect", xmin = 120, xmax = 180, ymin = -Inf, ymax = Inf,
           fill = "#C7E9C0", alpha = 0.35) +
  geom_vline(xintercept = 120, linetype = "dashed",
             color = "grey30", linewidth = 0.3) +
  geom_ridgeline(scale = 0.9, alpha = 0.35, linewidth = 0.4) +
  scale_fill_manual(values = lib_colors) +
  scale_color_manual(values = lib_colors) +
  scale_x_continuous(limits = c(35, 500), breaks = seq(50, 500, 50),
                     expand = c(0, 0)) +
  labs(x = "Fragment length (bp)", y = NULL,
       fill = "Library", color = "Library") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("analysis/frag_length_ridges.png", p_frag,
       width = 8, height = 0.35 * nrow(sample_order) + 1.2, dpi = 300)
message("Wrote analysis/frag_length_ridges.png")
