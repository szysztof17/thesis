#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette, make_sequential_cmap
import numpy as np
import re

# === Load style and config ===
load_custom_style()

palette = get_hue_palette(6)


method_order = ['PIDC', 'GENIE3', 'GRNBOOST2']
method_order = {'PIDC': 0,'GENIE3' : 1, 'GRNBOOST2' : 2}

cmap = make_sequential_cmap(sns.color_palette()[0])
#cmap = make_cmap_between(palette[0], palette[2])


# === File paths ===
base_name = "kim23_downsample_heatmap"
df = pd.read_csv(f'../data/{base_name}.csv', index_col=0)


dataset_names = df.columns.tolist()
pattern = re.compile(r"(\d+)pct-cells-(\d+)pct-reads")
cell_read_data = [pattern.search(name).groups() for name in dataset_names]



def plot_composite_figure(df, cell_read_data, dataset_name, base_name, cmap):
    """
    Plots a composite figure with a main heatmap and two marginal line plots.
    Layout: Heatmap takes 2/3 width, line plots take 1/3 width stacked vertically.
    """
    # === 1. Data Preparation ===
    values = df.loc[dataset_name].values
    plot_df = pd.DataFrame(cell_read_data, columns=["cells_pct", "reads_pct"])
    plot_df["cells_pct"] = plot_df["cells_pct"].astype(int)
    plot_df["reads_pct"] = plot_df["reads_pct"].astype(int)
    plot_df["AUROC"] = values

    # Data for the heatmap
    heatmap_data = plot_df.pivot(index="cells_pct", columns="reads_pct", values="AUROC").sort_index(ascending=False)
    
    # Data for the marginal plots
    marginal_reads_df = plot_df[plot_df['cells_pct'] == 100].sort_values('reads_pct')
    marginal_cells_df = plot_df[plot_df['reads_pct'] == 100].sort_values('cells_pct')

    # === 2. Create the Custom Figure Layout ===
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 3, figure=fig) # 2 rows, 3 columns grid

    # --- Heatmap (takes up 2/3 width and full height) ---
    ax1 = fig.add_subplot(gs[:, 0:2]) # Spans all rows, and columns 0 to 1
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        linewidths=0.5,
        cbar=True,
        cbar_kws={'label': 'AUROC', 'shrink': 0.8},
        ax=ax1,
        annot_kws={'size': 8, 'weight': 'bold', 'color': 'dimgray'}
    )
    cmap = ax1.collections[0].cmap
    norm = ax1.collections[0].norm

    for text in ax1.texts:
        value = float(text.get_text())
        # Map value to RGBA
        rgba = cmap(norm(value))
        # Compute perceived brightness (luma)
        brightness = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
        # Choose black or white depending on brightness
        text.set_color("#2e2d2c" if brightness > 0.5 else "white")
    ax1.set_title("AUROC Trade-off Matrix")
    ax1.set_xlabel(r"\% Reads")
    ax1.set_ylabel(r"\% Cells")
    
    # Get color limits to synchronize axes
    vmin, vmax = ax1.collections[0].get_clim()

    # --- Marginal Plot: Reads Effect (top-right) ---
    ax2 = fig.add_subplot(gs[0, 2]) # Row 0, column 2
    sns.lineplot(data=marginal_reads_df, x='reads_pct', y='AUROC', marker='o', ax=ax2)
    ax2.set_title(r"100\% Cells")
    ax2.set_xlabel(r"\% Reads")
    ax2.grid(True)
    ax2.set_ylim(vmin, vmax*1.02)
    ax2.set_xlim(-2,102)

    # --- Marginal Plot: Cells Effect (bottom-right) ---
    ax3 = fig.add_subplot(gs[1, 2]) # Row 1, column 2
    sns.lineplot(data=marginal_cells_df, x='cells_pct', y='AUROC', marker='o', ax=ax3)
    ax3.set_title(r"100\% Reads")
    ax3.set_xlabel(r"\% Cells")
    ax3.grid(True)
    ax3.set_ylim(vmin, vmax*1.02)
    ax3.set_xlim(-2,102)

    # === 3. Final Touches ===
    fig.suptitle(fr"{dataset_name} AUROC across \%Cells and \%Reads and marginal effects", fontsize='large')
    plt.tight_layout() # Adjust for main title

    plot_path = f"../plots/{base_name}_{dataset_name.lower()}_composite.pdf"
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()
    print(f"Saved figure: {plot_path}")


# === Run the plotting function for each dataset ===
for dataset in ["GENIE3", "PIDC"]:
    cmap = make_sequential_cmap(sns.color_palette()[method_order[dataset]])
    plot_composite_figure(df, cell_read_data, dataset, base_name, cmap)







