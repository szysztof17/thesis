#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette 
import numpy as np


# === Load style and config ===
load_custom_style()

palette = get_hue_palette(6)

# === File paths ===
base_name = "kim23_meta"
baseline_auroc = pd.read_csv(f'../data/kim23_meta_baseline.csv', index_col = 0).iloc[:, 0]  
long_df = pd.read_csv(f'../data/kim23_meta_data.csv', index_col = None)
plot_path = f"../plots/{base_name}.pdf"
print(long_df)

sns.lineplot(
    data=long_df,
    x="k", y="auroc",
    hue="method",
    style="metacells",
    markers=True,
    dashes=True,
    linewidth=1
)

# Extract the colors used by seaborn for each method
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
method_colors = {}
for handle, label in zip(handles, labels):
    if label in long_df['method'].unique():
        color = handle.get_color()
        method_colors[label] = color
print(method_colors)

for method, auroc_value in baseline_auroc.items():
    print(f'method: {method}, is in list: {method in method_colors}')
    if method in method_colors:
        print(f'loop for {method}')
        plt.axhline(
            y=auroc_value,
            color=method_colors[method],
            linestyle="dotted",
            linewidth=1.5,
            alpha=0.8
        )


plt.title("Effect of average cluster size (k) on AUROC")
plt.xlabel("Average cluster size (k) [log scale]")
plt.ylabel("AUROC")
plt.ylim((0.7,0.8))
plt.grid(True)
plt.xscale("log")
plt.legend()

baseline_line = mlines.Line2D(
    [], [], color='gray', linestyle='dotted', label='Baseline (no metacells)'
)
plt.legend(
    handles=[*plt.gca().get_legend_handles_labels()[0], baseline_line],
    loc="lower left" #, bbox_to_anchor=(1.05, 1)
)



set_figure_width(aspect_ratio = 2)
plt.tight_layout()

plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")




