"""
Plot ichorCNA summary metrics as horizontal bar charts.

Usage: python plot_tumor_fractions.py <tumor_fractions.csv> <samplesheet.csv> <output.pdf>
"""

import sys
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

csv_path, samplesheet_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

# ── load coverage from samplesheet ────────────────────────────────────────────
coverage_map = {}
with open(samplesheet_path) as f:
    for row in csv.DictReader(f):
        coverage_map[row["sample_id"]] = float(row["coverage"])

# ── load ichorCNA results ──────────────────────────────────────────────────────
samples, tf, ploidy, gc_mad, coverage, groups = [], [], [], [], [], []
with open(csv_path) as f:
    for row in csv.DictReader(f):
        s = row["sample"]
        samples.append(s)
        tf.append(float(row["tumor_frac"]))
        ploidy.append(float(row["ploidy"]))
        gc_mad.append(float(row["gc_mad"]))
        coverage.append(coverage_map.get(s, float("nan")))
        if s.startswith("IA"):
            groups.append("Autoimmune")
        elif s.startswith("IH") or s.startswith("BH"):
            groups.append("Healthy")
        else:
            groups.append("Cancer")

# sort by TF descending
order = sorted(range(len(tf)), key=lambda i: tf[i], reverse=True)
samples  = [samples[i]  for i in order]
tf       = [tf[i]       for i in order]
ploidy   = [ploidy[i]   for i in order]
gc_mad   = [gc_mad[i]   for i in order]
coverage = [coverage[i] for i in order]
groups   = [groups[i]   for i in order]

n = len(samples)
y = np.arange(n)

COLORS = {"Cancer": "#d7191c", "Autoimmune": "#2c7bb6", "Healthy": "#1a9641"}
bar_colors = [COLORS[g] for g in groups]

ROW_H = 0.18          # tighter rows
fig_h = n * ROW_H + 1.0
BAR_H = 0.72

fig, axes = plt.subplots(
    1, 4,
    figsize=(13, fig_h),
    sharey=True,
    gridspec_kw={"width_ratios": [2.2, 1.0, 1.0, 1.2]},
)
fig.subplots_adjust(wspace=0.07, left=0.11, right=0.97, top=0.95, bottom=0.05)

def style(ax):
    ax.grid(axis="x", which="major", color="#dddddd", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=7)
    ax.tick_params(axis="y", length=0)

# ── panel 1 : Tumor fraction (log scale) ──────────────────────────────────────
ax = axes[0]
ax.barh(y, tf, height=BAR_H, color=bar_colors, linewidth=0)
ax.set_xscale("log")
ax.set_xlabel("Tumor fraction", fontsize=8)
ax.set_title("Tumor\nfraction", fontsize=9, fontweight="bold", linespacing=1.2)
ax.set_yticks(y)
ax.set_yticklabels(samples, fontsize=6)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
style(ax)

# ── panel 2 : Ploidy ──────────────────────────────────────────────────────────
ax = axes[1]
ax.barh(y, ploidy, height=BAR_H, color=bar_colors, linewidth=0)
ax.axvline(2.0, color="#444444", linewidth=0.8, linestyle="--", zorder=3)
ax.set_xlabel("Ploidy", fontsize=8)
ax.set_title("Ploidy", fontsize=9, fontweight="bold")
pad = 0.15
ax.set_xlim(min(ploidy) - pad, max(ploidy) + pad)
style(ax)

# ── panel 3 : GC-map MAD ──────────────────────────────────────────────────────
ax = axes[2]
ax.barh(y, gc_mad, height=BAR_H, color=bar_colors, linewidth=0)
ax.axvline(0.05, color="#e08020", linewidth=0.8, linestyle="--", zorder=3)
ax.set_xlabel("GC-map MAD", fontsize=8)
ax.set_title("GC-map\nMAD", fontsize=9, fontweight="bold", linespacing=1.2)
ax.set_xlim(0, max(gc_mad) * 1.15)
style(ax)

# ── panel 4 : Coverage ────────────────────────────────────────────────────────
ax = axes[3]
ax.barh(y, coverage, height=BAR_H, color=bar_colors, linewidth=0)
ax.set_xscale("log")
ax.set_xlabel("Coverage (×)", fontsize=8)
ax.set_title("Coverage", fontsize=9, fontweight="bold")
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
style(ax)

# ── legend ────────────────────────────────────────────────────────────────────
legend_patches = [mpatches.Patch(color=c, label=g) for g, c in COLORS.items()]
fig.legend(handles=legend_patches, loc="lower center", ncol=3,
           fontsize=8, framealpha=0.8, bbox_to_anchor=(0.54, 0.0))

for ax in axes:
    ax.invert_yaxis()

plt.savefig(out_path, bbox_inches="tight", dpi=150)
print(f"Saved {out_path}")
