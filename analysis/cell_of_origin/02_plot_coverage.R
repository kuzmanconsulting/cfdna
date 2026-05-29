#!/usr/bin/env Rscript
# Cell-of-origin cfDNA coverage wiggle.
# Reads the per-(sample,tissue) offset profiles from 01_coverage.sh, read-scales each profile by
# library size, and overlays each cancer sample on its matched normal at every tissue's
# super-enhancers. A deeper central dip (offset 0) in the cancer than in its matched normal at a
# tissue's SEs means that tissue sheds more cfDNA -> a cell-of-origin signal. Run from repo root.

suppressPackageStartupMessages(library(tidyverse))

cov_dir <- "analysis/cell_of_origin/coverage"
out_dir <- "analysis/cell_of_origin"
group_colors <- c(healthy = "#523FCC", cancer = "#FF004C")

# ---- load offset profiles; mean cov per offset, read-scaled by library size -------------
read_one <- function(f) {
  h <- readLines(f, n = 1)
  n_sites <- as.integer(str_match(h, "n_sites=(\\d+)")[, 2])
  reads   <- as.numeric(str_match(h, "total_reads=(\\d+)")[, 2])
  bn <- str_match(basename(f), "^(.*)\\.([^.]+)\\.offset_cov\\.tsv$")
  read_tsv(f, comment = "#", show_col_types = FALSE) |>
    transmute(sample = bn[, 2], tissue = bn[, 3], offset,
              scaled = (sum_cov / n_sites) / reads * 1e6)   # mean cov per Mread
}
files <- list.files(cov_dir, pattern = "\\.offset_cov\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No profiles in ", cov_dir, " -- run 01_coverage.sh first.")
prof <- map(files, read_one) |> list_rbind()

# ---- metadata: group, own tissue, matched-normal stratum --------------------------------
own <- read_csv("analysis/tumor_types.csv", show_col_types = FALSE) |>
  select(disease, own_tissue = tissue)
meta <- read_csv("snyder2016_metadata_GSE.csv", show_col_types = FALSE) |>
  select(sample_id, sample_group, disease, library_type, coverage) |>
  left_join(own, by = "disease") |>
  mutate(stratum = case_when(library_type == "DSP" ~ "DSP",
                             coverage >= 10        ~ "SSP_ge10",
                             TRUE                  ~ "SSP_lt10"))
ctrl_of <- c(DSP = "IH01", SSP_ge10 = "IH02", SSP_lt10 = "IH03")
prof <- prof |> left_join(meta, by = c("sample" = "sample_id"))

# ---- cancer vs matched normal: overlaid scaled-coverage wiggle, faceted -----------------
normal <- prof |> filter(sample %in% ctrl_of) |>
  select(ctrl = sample, tissue, offset, n_scaled = scaled)

panels <- prof |> filter(sample_group == "cancer") |>
  mutate(ctrl = ctrl_of[stratum]) |>
  left_join(normal, by = c("ctrl", "tissue", "offset")) |>
  pivot_longer(c(scaled, n_scaled), names_to = "role", values_to = "cov") |>
  mutate(role = if_else(role == "scaled", "cancer", "matched normal"))

for (s in unique(panels$sample)) {
  pd <- filter(panels, sample == s)
  ot <- pd$own_tissue[1]
  p <- ggplot(pd, aes(offset, cov, colour = role)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ tissue, scales = "free_y") +
    scale_colour_manual(values = c(cancer = "#FF004C", `matched normal` = "#523FCC"), name = NULL) +
    labs(title = paste0(s, " (own tissue: ", ot, ")"),
         x = "Position relative to super-enhancer center (bp)",
         y = "Coverage per million reads") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(out_dir, paste0("wiggle_", s, ".png")), p, width = 10, height = 6, dpi = 300)
  message("Wrote wiggle_", s, ".png")
}

# ---- COO summary: cancer/normal central-dip ratio per tissue ----------------------------
# dip = 1 - mean(cov[|offset|<=150]) / mean(cov[|offset|>=750]); ratio = cancer_dip / normal_dip
dip <- prof |> group_by(sample, tissue, sample_group, stratum, own_tissue) |>
  summarise(dip = 1 - mean(scaled[abs(offset) <= 150]) / mean(scaled[abs(offset) >= 750]),
            .groups = "drop")
ndip <- dip |> filter(sample %in% ctrl_of) |> select(ctrl = sample, tissue, n_dip = dip)
score <- dip |> filter(sample_group == "cancer") |>
  mutate(ctrl = ctrl_of[stratum]) |>
  left_join(ndip, by = c("ctrl", "tissue")) |>
  mutate(dip_ratio = dip / n_dip, is_own = !is.na(own_tissue) & own_tissue == tissue)
write_tsv(score, file.path(out_dir, "coo_dip_ratio.tsv"))

p_hm <- ggplot(score, aes(tissue, sample, fill = dip_ratio)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = filter(score, is_own), colour = "black", linewidth = 0.8, fill = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 1, name = "dip ratio\n(cancer/normal)") +
  labs(title = "Super-enhancer central-dip ratio (black = sample's own tissue)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "coo_dip_heatmap.png"), p_hm, width = 7, height = 6, dpi = 300)
message("Wrote coo_dip_heatmap.png; done.")
