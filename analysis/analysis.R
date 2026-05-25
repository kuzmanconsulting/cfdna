#!/usr/bin/env Rscript
# Load samples from cfdna-finale-snakemake/samplesheet.csv, attach disease
# (from snyder2016_metadata_GSE.csv) and disease_label (from
# analysis/tumor_types.csv).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
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
