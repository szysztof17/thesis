#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette 

# === Load style and config ===
load_custom_style()

# === File paths ===
base_name = "QC_influence"
df_counts = pd.read_csv('../data/QC_nGenesCounts_filtering_new.csv')

df_percentile1 = pd.read_csv('../data/QC_percentile_based_filtering.csv')
plot_path = f"../plots/{base_name}.pdf"

df_percentile1["filter_type"] = df_percentile1["percentile_high"].apply(
    lambda x: "low-high" if x < 100 else "low-100"
)
df_percentile = df_percentile1[df_percentile1["filter_type"] != "low-100"]
df_percentile2 = df_percentile1[df_percentile1["mt"] == 25]
# === Create consistent hue palette ===
all_mt_values = sorted(set(df_counts["mt"]).union(df_percentile["mt"]), reverse=True)
colors = get_hue_palette(len(all_mt_values))
colors[0] = '#BAA898'
colors[2] = '#73ba9b'
palette = dict(zip(all_mt_values, colors))

palette2 = get_hue_palette(len(all_mt_values) + len(set(df_percentile1["filter_type"])))
palette2 = palette2[len(all_mt_values):]
fig, axes = plt.subplots(1, 3)  

# === PLOT 1: Percentile-based filtering ===

sns.lineplot(
    data=df_counts.sort_values("n_genes"),
    x="n_genes",
    y="AUROC",
    hue="mt",
    ax=axes[0],
    marker="o",
    markersize=4,

    palette=palette
)
axes[0].set_title("nGenesCounts filtering")
axes[0].set_xlabel("Minimum n_genes_by_counts")
axes[0].set_ylabel("AUROC")
legend = axes[0].legend(title=r"MT cutoff (\%)", fontsize='x-small')
legend.get_title().set_fontsize('small')

# === PLOT 2: Fixed nGenesCounts filtering ===
sns.lineplot(
    data=df_percentile,
    x="percentile_low",
    y="AUROC",
    hue="mt",
    markersize=4,

    marker="o",
    ax=axes[1],
    palette=palette,
)
axes[1].set_title("Two-sided pct filtering")
axes[1].set_xlabel("Lower and upper pct threshold")
axes[1].set_ylabel("AUROC")
axes[1].legend(title=r"MT cutoff (\%)", fontsize='x-small')
legend = axes[1].legend(title=r"MT cutoff (\%)", fontsize='x-small')
legend.get_title().set_fontsize('small')

axes[0].set_ylim(0.5, 1.0)
axes[1].set_ylim(0.5, 1.0)

axes[0].axvline(x=1259, color='red', linestyle='--', linewidth=1, alpha=0.5)
axes[0].text(
    1259,                             # x-position (same as the line)
    axes[1].get_ylim()[1] * 0.9,            # y-position at the bottom of the y-axis
    "1259: local minimum",
    color="gray",
    fontsize=7,
    ha="left",                       # Horizontal alignment: left of the x-position
    va="bottom",                      # Vertical alignment: bottom of the y-position
    alpha=0.7
)
# Set common y-axis limits

# === PLOT 3: Compare low-high vs low-100 ===

df_percentile2.loc[:, "filter_type"] = df_percentile2["filter_type"].replace({
    "low-high": "Two-sided",
    "low-100": "One-sided"
})

sns.lineplot(
    data=df_percentile2,
    x="percentile_low",
    y="AUROC",
    hue="filter_type",
    markersize=4,

    marker="o",
    ax=axes[2],
    palette = palette2
)
axes[2].set_title("One vs two-sided pct filtering")
axes[2].set_xlabel("Lower pct threshold")
axes[2].set_ylabel("AUROC")
axes[2].set_ylim(0.5, 1.0)
legend = axes[2].legend(title="Filter type")
legend.get_title().set_fontsize('small')

axes[1].set_ylabel("")
axes[1].tick_params(labelleft=False)

axes[2].set_ylabel("")
axes[2].tick_params(labelleft=False)

# Adjust layout
#plt.subplots_adjust(wspace=0.25)

set_figure_width(aspect_ratio=2.25)
plt.tight_layout()

plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
