"""Plot fragment length distribution, annotated with the analytical length classes.

QC/context figure built from `md.bam` via `finaletoolkit frag-length-bins
--bin-size 1` (per-length counts). It is independent of `length_classes` in
config.yaml; the shaded bands are annotations of classes III (S-WPS 35-80)
and IV (L-WPS 120-180), and the dashed line marks the I | II split at 120.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

tsv_path = snakemake.input.tsv
png_path = snakemake.output.png
sample = snakemake.wildcards.sample

df = pd.read_csv(tsv_path, sep="\t", comment="#")

# FinaleToolkit frag-length-bins schema: min, max, count (one row per bin).
# Fragment length = bin midpoint (== min == max when --bin-size 1).
df = df.sort_values("min")
lengths = ((df["min"] + df["max"]) / 2).to_numpy()
counts = df["count"].to_numpy()

mask_III = (lengths >= 35) & (lengths <= 80)    # Snyder strict S-WPS
mask_IV = (lengths >= 120) & (lengths <= 180)    # Snyder strict L-WPS

fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=False)

ax = axes[0]
ax.fill_between(lengths, counts, step="mid", alpha=0.7)
ax.axvspan(35, 80, color="tab:orange", alpha=0.15, label="Class III (S-WPS 35-80)")
ax.axvspan(120, 180, color="tab:green", alpha=0.15, label="Class IV (L-WPS 120-180)")
ax.axvline(120, color="tab:red", linestyle="--", linewidth=1.2, label="I | II boundary (120)")
ax.set_xlim(0, max(400, lengths.max()))
ax.set_xlabel("Fragment length (bp)")
ax.set_ylabel("Fragment count")
ax.set_title(f"{sample} — full fragment-length distribution")
ax.legend(loc="upper right")

ax = axes[1]
ax.bar(lengths[mask_III], counts[mask_III], width=1.0, color="tab:orange")
ax.set_xlabel("Fragment length (bp)")
ax.set_ylabel("Fragment count")
ax.set_title(f"{sample} — Class III (S-WPS, 35-80 bp)")

ax = axes[2]
ax.bar(lengths[mask_IV], counts[mask_IV], width=1.0, color="tab:green")
ax.set_xlabel("Fragment length (bp)")
ax.set_ylabel("Fragment count")
ax.set_title(f"{sample} — Class IV (L-WPS, 120-180 bp)")

fig.tight_layout()
fig.savefig(png_path, dpi=120)
