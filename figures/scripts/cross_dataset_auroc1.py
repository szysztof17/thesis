#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette 
import numpy as np


# === Load style and config ===
load_custom_style()

palette = get_hue_palette(6, style = 'two')

# === File paths ===
base_name = "cross_dataset_auroc"
df = pd.read_csv(f'../data/{base_name}.csv')
plot_path = f"../plots/{base_name}1.pdf"

# Plot configuration

entries = {
    "hhep": {
        "500": "hhep_STRING12_500TFs",
        "1000": "hhep_STRING12_1000TFs"
    },
    "kim23-hm1-extgenes": {
        "500": "kim23-hm1-extgenes500",
        "1000": "kim23-hm1-extgenes1000"
    }
}

legend_labels = {
    "hhep": "hHep",
    "kim23-hm1-extgenes": "Kim23-hm1-extgenes" 
}

# Prepare data
sizes = ["500", "1000"]
x = np.arange(len(sizes))
width = 0.35

fig, ax = plt.subplots()

for i, dataset in enumerate(entries.keys()):
    medians = []
    err_low = []
    err_high = []

    for size in sizes:
        col = entries[dataset][size]
        values = df[col]
        medians.append(values.median())
        err_low.append(values.median() - values.min())
        err_high.append(values.max() - values.median())
    
    # Plot bars
    bar_positions = x + (i - 0.5) * width
    bars = ax.bar(bar_positions, medians, width=width,
                  yerr=[err_low, err_high], capsize=6,
                  label=legend_labels.get(dataset, dataset), color=palette[i])

    # Add value labels on top of each bar
    for xpos, value in zip(bar_positions, medians):
        ax.text(xpos, value + 0.03, f"{value:.3f}",
                ha='center', va='bottom', fontsize=9)

ax.set_xticks(x)
ax.set_xticklabels(sizes)
ax.set_ylabel("AUROC")
ax.set_xlabel("Gene set size")
ax.set_ylim(0.5, 1)
ax.set_title("AUROC on shared gene set (hHep-based)")
ax.legend(title="Dataset")

set_figure_width(aspect_ratio=1, n_per_page=2)
plt.tight_layout()

plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")







# # Prepare data
# labels = list(entries.keys())
# x = np.arange(len(labels))
# width = 0.35

# fig, ax = plt.subplots(figsize=(8, 5))

# for i, size in enumerate(["500", "1000"]):
#     medians = []
#     err_low = []
#     err_high = []

#     for dataset in labels:
#         col = entries[dataset][size]
#         values = df[col]
#         medians.append(values.median())
#         err_low.append(values.median() - values.min())
#         err_high.append(values.max() - values.median())

#     ax.bar(x + (i - 0.5)*width, medians, width=width,
#            yerr=[err_low, err_high], capsize=6, label=f"{size} TFs", color=palette[i+3])

# ax.set_xticks(x)
# ax.set_xticklabels(labels)
# ax.set_ylabel("AUROC")
# ax.set_ylim(0.5, 1)
# ax.set_title("AUROC on shared gene set (hhep-based)")
# ax.legend(title="Gene set size")




# set_figure_width(aspect_ratio=1, n_per_page = 2)
# plt.tight_layout()

# plt.savefig(plot_path, bbox_inches='tight')
# print(f"Save a figure to: {plot_path}")
