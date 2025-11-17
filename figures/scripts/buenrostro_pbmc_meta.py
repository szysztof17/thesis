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


base_name = 'buenrostro_pbmc_meta'
# 1. Load Data
pbmc_df = pd.read_csv('../data/pbmc10k_meta-auroc.csv', index_col=0)
buen_df = pd.read_csv('../data/buenrostro18-auroc.csv', index_col=0)

plot_path = f"../plots/{base_name}.pdf"


def process_df(df):


    if 'pbmc10k_meta_STRING12_1000TFs' in df.columns:
        baseline_col = 'pbmc10k_meta_STRING12_1000TFs'
    elif 'buenrostro18_STRING12_1000TFs' in df.columns:
        baseline_col = 'buenrostro18_STRING12_1000TFs'
    else:
        # 2. fallback: find the first column that doesn't contain aggregation keywords
        # This is robust for future datasets
        try:
            baseline_col = [c for c in df.columns 
                            if 'KMeans' not in c 
                            and 'SEACells' not in c 
                            and 'pct' not in c][0]
        except IndexError:
            raise ValueError("Could not identify a baseline column in the dataframe.")


    baseline_series = df[baseline_col]
    
    # Filter & Melt Metacell columns
    metacell_cols = [c for c in df.columns if 'KMeans' in c or 'SEACells' in c]
    df_meta = df[metacell_cols].copy().reset_index().rename(columns={'index': 'grn_method'})
    df_long = df_meta.melt(id_vars='grn_method', var_name='sample_name', value_name='AUROC')
    
    # Extract N_metacells (nmc) and Average_size (avg)
    regex = r'(?P<method>KMeans|SEACells)-nmc(?P<nmc>\d+)-avg(?P<avg>\d+)'
    extracted = df_long['sample_name'].str.extract(regex)
    df_long = pd.concat([df_long, extracted], axis=1)
    df_long[['nmc', 'avg']] = df_long[['nmc', 'avg']].apply(pd.to_numeric)
    
    return df_long, baseline_series

pbmc_long, pbmc_base = process_df(pbmc_df)
buen_long, buen_base = process_df(buen_df)

# 3. Plotting Side-by-Side
fig, axes = plt.subplots(1, 2, sharey=True)

def plot_panel(ax, df, baseline, title):
    sns.lineplot(
        data=df, x="avg", y="AUROC", hue="grn_method", style="method",
        markers=True, dashes=True, linewidth=1, palette=palette, ax=ax
    )
    # Add Baselines
    for method, score in baseline.items():
        if method in palette:
            ax.axhline(y=score, color=palette[method], linestyle='dotted', linewidth=1.5, alpha=0.8)
            
    ax.set_title(title)
    ax.set_xscale("log")
    ax.grid(True, which='major', linestyle='--', alpha=0.5)
    ax.get_legend().remove() # Clean up individual legends
    ax.set_xlabel("Average cluster size (k) [log scale]")

plot_panel(axes[1], pbmc_long, pbmc_base, "PBMC10k ")
plot_panel(axes[0], buen_long, buen_base, "Buenrostro18 ")

axes[0].set_ylabel("AUROC")

# 4. Unified Legend
handles, labels = axes[0].get_legend_handles_labels()
handles.append(mlines.Line2D([], [], color='black', linestyle=':', label='Baseline'))
labels.append('Baseline')
axes[0].legend(
    handles,
    labels,
    loc='lower left',
    #fontsize=8,                   # small text
    fontsize='x-small',                   # small text    
    labelspacing=0.3,             # tighten vertical spacing
    handlelength=1.2,             # shorter line samples
    handletextpad=0.4             # reduce space between line and text
)
set_figure_width(aspect_ratio = 2)
plt.tight_layout()   # reserve right margin for legend

#plt.subplots_adjust(bottom=0.15) # Space for legend
#set_figure_width(aspect_ratio = 2)

plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
