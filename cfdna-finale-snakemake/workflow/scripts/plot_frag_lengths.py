"""Plot fragment length distribution with overlays for S-WPS (35-80) and L-WPS (120-180).

Consumes a TSV produced by `finaletoolkit frag-length-bins --bin-size 1`,
which contains per-length fragment counts.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

tsv_path = snakemake.input.tsv
png_path = snakemake.output.png
sample = snakemake.wildcards.sample

df = pd.read_csv(tsv_path, sep="\t", comment="#")

# Heuristic column detection: first numeric column is length, second is count.
length_col = next(c for c in df.columns if "length" in c.lower() or "len" in c.lower())
count_col = next(c for c in df.columns if c != length_col and pd.api.types.is_numeric_dtype(df[c]))

df = df.sort_values(length_col)
lengths = df[length_col].to_numpy()
counts = df[count_col].to_numpy()

mask_short = (lengths >= 35) & (lengths <= 80)
mask_long = (lengths >= 120) & (lengths <= 180)

fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=False)

ax = axes[0]
ax.fill_between(lengths, counts, step="mid", alpha=0.7)
ax.axvspan(35, 80, color="tab:orange", alpha=0.15, label="S-WPS 35-80")
ax.axvspan(120, 180, color="tab:green", alpha=0.15, label="L-WPS 120-180")
ax.set_xlim(0, max(400, lengths.max()))
ax.set_xlabel("Fragment length (bp)")
ax.set_ylabel("Fragment count")
ax.set_title(f"{sample} — full fragment-length distribution")
ax.legend(loc="upper right")

ax = axes[1]
ax.bar(lengths[mask_short], counts[mask_short], width=1.0, color="tab:orange")
ax.set_xlabel("Fragment length (bp)")
ax.set_ylabel("Fragment count")
ax.set_title(f"{sample} — short fragments (35-80 bp, S-WPS)")

ax = axes[2]
ax.bar(lengths[mask_long], counts[mask_long], width=1.0, color="tab:green")
ax.set_xlabel("Fragment length (bp)")
ax.set_ylabel("Fragment count")
ax.set_title(f"{sample} — long fragments (120-180 bp, L-WPS)")

fig.tight_layout()
fig.savefig(png_path, dpi=120)
