#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette, METHOD_PALETTE_INDICES, get_method_color

import numpy as np
import re

# === Load style and config ===
load_custom_style()

#palette = get_hue_palette(6)
palette = {m: get_method_color(m) for m in METHOD_PALETTE_INDICES}

# === File paths ===
base_name = "pbmc10k_meta"
pbmc10k_auroc_file = '../data/pbmc10k_meta-auroc.csv'
pbmc10k_auroc = pd.read_csv(pbmc10k_auroc_file, index_col=0)
plot_path = f"../plots/{base_name}.pdf"

df = pbmc10k_auroc.copy()
df_basic = df[['pbmc10k_meta_STRING12_1000TFs']]

metacells_cols = df.columns[df.columns.str.contains('KMeans|SEACells')]
df_metacells = df[metacells_cols]


df_metacells.index.name = 'grn_method' 
df_metacells = df_metacells.reset_index()
df_long = df_metacells.melt(var_name='sample_name', value_name='AUROC', id_vars=['grn_method'])
regex_pattern = r'pbmc10k_meta-(?P<method>KMeans|SEACells)-nmc(?P<nmc>\d+)-avg(?P<avg>\d+)'
df_long[['method', 'nmc', 'avg']] = df_long['sample_name'].str.extract(regex_pattern)
df_long['nmc'] = pd.to_numeric(df_long['nmc'])
df_long['avg'] = pd.to_numeric(df_long['avg'])
df_long = df_long.drop(columns=['sample_name'])


baseline_auroc = df_basic.iloc[:, 0]  
long_df = df_long





sns.lineplot(
    data=df_long,
    x="avg", y="AUROC",
    hue="grn_method",
    style="method",
    markers=True,
    dashes=True,
    linewidth=1,
    palette=palette,
)


for method, auroc_value in baseline_auroc.items():
    print(f'method: {method}, is in list: {method in palette}')
    if method in palette:
        print(f'loop for {method}')
        plt.axhline(
            y=auroc_value,
            color=palette[method],
            linestyle="dotted",
            linewidth=1.5,
            alpha=0.8
        )


plt.title("Effect of average cluster size (k) on AUROC [pbmc10k data]")
plt.xlabel("Average cluster size (k) [log scale]")
plt.ylabel("AUROC")
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




