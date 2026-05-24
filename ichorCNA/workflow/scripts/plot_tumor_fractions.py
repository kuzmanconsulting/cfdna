"""
Plot ichorCNA summary metrics as horizontal bar charts.

Usage: python plot_tumor_fractions.py <tumor_fractions.csv> <output.pdf>
"""

import sys
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

csv_path, out_path = sys.argv[1], sys.argv[2]

# ── load ──────────────────────────────────────────────────────────────────────
samples, tf, ploidy, gc_mad, groups = [], [], [], [], []
with open(csv_path) as f:
    for row in csv.DictReader(f):
        samples.append(row["sample"])
        tf.append(float(row["tumor_frac"]))
        ploidy.append(float(row["ploidy"]))
        gc_mad.append(float(row["gc_mad"]))
        groups.append("Autoimmune" if row["sample"].startswith("IA") else "Cancer")

# sort by TF descending
order = sorted(range(len(tf)), key=lambda i: tf[i], reverse=True)
samples = [samples[i] for i in order]
tf      = [tf[i]      for i in order]
ploidy  = [ploidy[i]  for i in order]
gc_mad  = [gc_mad[i]  for i in order]
groups  = [groups[i]  for i in order]

n = len(samples)
y = np.arange(n)

COLORS = {"Cancer": "#d7191c", "Autoimmune": "#2c7bb6"}
bar_colors = [COLORS[g] for g in groups]

ROW_H = 0.22          # inches per sample row
fig_h = n * ROW_H + 1.2   # +margin for title/labels
BAR_H = 0.7           # bar height in data units (out of 1.0)

fig, axes = plt.subplots(
    1, 3,
    figsize=(12, fig_h),
    sharey=True,
    gridspec_kw={"width_ratios": [2.5, 1.2, 1.2]},
)
fig.subplots_adjust(wspace=0.08, left=0.13, right=0.97, top=0.96, bottom=0.06)

# ── panel 1 : Tumor fraction (log scale) ──────────────────────────────────────
ax = axes[0]
ax.barh(y, tf, height=BAR_H, color=bar_colors, linewidth=0)
ax.set_xscale("log")
ax.set_xlabel("Tumor fraction (log scale)", fontsize=9)
ax.set_title("Tumor fraction", fontsize=10, fontweight="bold")
ax.set_yticks(y)
ax.set_yticklabels(samples, fontsize=6.5)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.tick_params(axis="x", labelsize=8)
ax.grid(axis="x", which="major", color="#dddddd", linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

# ── panel 2 : Ploidy ──────────────────────────────────────────────────────────
ax = axes[1]
ax.barh(y, ploidy, height=BAR_H, color=bar_colors, linewidth=0)
ax.axvline(2.0, color="#444444", linewidth=0.9, linestyle="--", zorder=3,
           label="diploid (2.0)")
ax.set_xlabel("Ploidy", fontsize=9)
ax.set_title("Ploidy", fontsize=10, fontweight="bold")
ax.tick_params(axis="x", labelsize=8)
ax.grid(axis="x", which="major", color="#dddddd", linewidth=0.5, zorder=0)
ax.set_axisbelow(True)
# narrow x range to show variation
pad = 0.15
ax.set_xlim(min(ploidy) - pad, max(ploidy) + pad)
ax.legend(fontsize=7, loc="lower right", framealpha=0.7)

# ── panel 3 : GC-map MAD ──────────────────────────────────────────────────────
ax = axes[2]
ax.barh(y, gc_mad, height=BAR_H, color=bar_colors, linewidth=0)
ax.axvline(0.05, color="#e08020", linewidth=0.9, linestyle="--", zorder=3,
           label="QC threshold (0.05)")
ax.set_xlabel("GC-map MAD", fontsize=9)
ax.set_title("GC-map MAD", fontsize=10, fontweight="bold")
ax.tick_params(axis="x", labelsize=8)
ax.grid(axis="x", which="major", color="#dddddd", linewidth=0.5, zorder=0)
ax.set_axisbelow(True)
ax.set_xlim(0, max(gc_mad) * 1.15)
ax.legend(fontsize=7, loc="lower right", framealpha=0.7)

# ── shared legend ─────────────────────────────────────────────────────────────
legend_patches = [mpatches.Patch(color=c, label=g) for g, c in COLORS.items()]
fig.legend(handles=legend_patches, loc="lower center", ncol=2,
           fontsize=8, framealpha=0.8,
           bbox_to_anchor=(0.55, 0.005))

# invert Y so highest-TF sample is at top
for ax in axes:
    ax.invert_yaxis()

plt.savefig(out_path, bbox_inches="tight", dpi=150)
print(f"Saved {out_path}")
