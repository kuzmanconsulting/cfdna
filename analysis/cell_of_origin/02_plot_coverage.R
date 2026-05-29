#!/usr/bin/env Rscript
# Cell-of-origin cfDNA coverage analysis.
# Reads the per-(sample,tissue) offset profiles from 01_coverage.sh, flank-normalizes,
# and compares each cancer sample to its matched-by-stratum healthy control at
# tissue-specific TFBS sets. A deeper central coverage dip at a tissue's regulatory sites
# implies that tissue contributes more cfDNA -> a tissue-of-origin signal.
#
# Run from the repo root (uses relative paths).

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

cov_dir  <- "analysis/cell_of_origin/coverage"
out_dir  <- "analysis/cell_of_origin"
cache    <- "analysis/.cache/coo_offset_cov.tsv.gz"
dir.create("analysis/.cache", showWarnings = FALSE, recursive = TRUE)

group_colors <- c(healthy = "#523FCC", cancer = "#FF004C")

# sample -> tissue of origin (from Snyder clinical dx); healthy samples have none
sample_tissue <- tribble(
  ~sample_id, ~own_tissue,
  "IC37","colon", "IC33","colon",
  "IC15","lung",  "IC32","lung", "IC20","lung", "IC28","lung", "IC10","lung",
  "IC17","liver", "IC23","liver",
  "IC35","breast","IC46","breast",
  "IC49","pancreas","IC50","pancreas","IC51","pancreas","IC52","pancreas"
)

# ---- load offset profiles (cache to avoid re-reading) -----------------------------------
read_one <- function(f) {
  n_sites <- as.integer(str_match(readLines(f, n = 1), "n_sites=(\\d+)")[, 2])
  bn <- str_match(basename(f), "^(.*)\\.([^.]+)\\.offset_cov\\.tsv$")
  read_tsv(f, comment = "#", show_col_types = FALSE) |>
    transmute(sample_id = bn[, 2], tissue = bn[, 3],
              offset, mean_cov = sum_cov / n_sites, n_sites = n_sites)
}

if (file.exists(cache)) {
  prof <- read_tsv(cache, show_col_types = FALSE)
} else {
  files <- list.files(cov_dir, pattern = "\\.offset_cov\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("No coverage profiles in ", cov_dir, " — run 01_coverage.sh first.")
  prof <- map(files, read_one) |> list_rbind()
  write_tsv(prof, cache)
}
message("Loaded ", n_distinct(prof$sample_id), " samples x ", n_distinct(prof$tissue), " tissues")

# ---- sample metadata + matched-control stratum ------------------------------------------
samples <- read_csv("cfdna-finale-snakemake/samplesheet.csv", show_col_types = FALSE) |>
  mutate(stratum = case_when(
    library_type == "DSP"                  ~ "DSP",
    library_type == "SSP" & coverage >= 10 ~ "SSP_ge10",
    TRUE                                   ~ "SSP_lt10"))
stratum_ctrl <- c(DSP = "IH01", SSP_ge10 = "IH02", SSP_lt10 = "IH03")

prof <- prof |>
  left_join(select(samples, sample_id, sample_group, library_type, coverage, stratum),
            by = "sample_id") |>
  left_join(sample_tissue, by = "sample_id")

# ---- flank normalization: relcov = mean_cov / mean(mean_cov over |offset| in [750,1000]) -
prof <- prof |>
  group_by(sample_id, tissue) |>
  mutate(flank = mean(mean_cov[abs(offset) >= 750]),
         relcov = mean_cov / flank) |>
  ungroup()

# ---- per-tissue profiles: each cancer sample vs its matched healthy control --------------
cancers <- samples |> filter(sample_group == "cancer")
tissues <- sort(unique(prof$tissue))

for (ti in tissues) {
  pd <- map(cancers$sample_id, function(s) {
    ctrl <- stratum_ctrl[[ samples$stratum[samples$sample_id == s] ]]
    prof |>
      filter(tissue == ti, sample_id %in% c(s, ctrl)) |>
      mutate(panel = s,
             role  = if_else(sample_id == s, "cancer", "matched normal"))
  }) |> list_rbind()
  if (nrow(pd) == 0) next

  p <- ggplot(pd, aes(offset, relcov, color = sample_group, group = sample_id)) +
    geom_hline(yintercept = 1, linewidth = 0.25, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_line(linewidth = 0.5) +
    facet_wrap(~panel) +
    scale_color_manual(values = group_colors, name = NULL) +
    labs(title = paste0("Coverage at ", ti, " TFBS (cancer vs matched normal)"),
         x = "Position relative to TFBS center (bp)", y = "Normalized coverage") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(out_dir, paste0("profiles_", ti, ".png")), p,
         width = 11, height = 7, dpi = 300)
  message("Wrote profiles_", ti, ".png")
}

# ---- dip depth = 1 - mean(relcov, |offset| <= 150) ---------------------------------------
dip <- prof |>
  group_by(sample_id, tissue, sample_group, stratum, own_tissue) |>
  summarise(dip = 1 - mean(relcov[abs(offset) <= 150]), .groups = "drop")

# heatmap: sample x tissue dip depth; own-tissue cells outlined
samp_order <- samples |>
  left_join(sample_tissue, by = "sample_id") |>
  arrange(sample_group, own_tissue, sample_id) |> pull(sample_id)
dip <- dip |> mutate(sample_id = factor(sample_id, levels = samp_order),
                     is_own = !is.na(own_tissue) & own_tissue == tissue)

p_hm <- ggplot(dip, aes(tissue, sample_id, fill = dip)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(data = filter(dip, is_own), color = "black", linewidth = 0.8, fill = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "dip depth") +
  labs(title = "Coverage dip depth by sample x tissue (black = sample's own tissue)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "coo_dipdepth_heatmap.png"), p_hm,
       width = 7, height = 8, dpi = 300)
message("Wrote coo_dipdepth_heatmap.png")

# delta dip at each cancer sample's OWN tissue: cancer minus its matched normal
delta <- cancers |>
  left_join(sample_tissue, by = "sample_id") |>
  filter(!is.na(own_tissue)) |>
  rowwise() |>
  mutate(ctrl = stratum_ctrl[[stratum]],
         dip_cancer = dip$dip[dip$sample_id == sample_id & dip$tissue == own_tissue][1],
         dip_normal = dip$dip[dip$sample_id == ctrl      & dip$tissue == own_tissue][1],
         delta_dip  = dip_cancer - dip_normal) |>
  ungroup()

p_delta <- ggplot(delta, aes(reorder(sample_id, delta_dip), delta_dip, fill = own_tissue)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  coord_flip() +
  labs(title = "Cancer minus matched-normal dip at the sample's own tissue",
       x = NULL, y = "delta dip depth (cancer - normal)", fill = "tissue") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "coo_dip_delta.png"), p_delta, width = 7, height = 5, dpi = 300)
message("Wrote coo_dip_delta.png")

write_tsv(dip, file.path(out_dir, "dip_depth.tsv"))
message("Done.")
