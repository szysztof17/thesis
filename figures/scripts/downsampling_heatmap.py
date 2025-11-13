#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette, make_sequential_cmap
import numpy as np
import re

# === Load style and config ===
load_custom_style()

palette = get_hue_palette(6)

cmap = make_sequential_cmap(sns.color_palette()[1])
#cmap = make_cmap_between(palette[0], palette[2])


# === File paths ===
base_name = "kim23_downsample_heatmap"
df = pd.read_csv(f'../data/{base_name}.csv', index_col=0)


dataset_names = df.columns.tolist()
pattern = re.compile(r"(\d+)pct-cells-(\d+)pct-reads")
cell_read_data = [pattern.search(name).groups() for name in dataset_names]



def plot_auroc_heatmap(df, cell_read_data, dataset_name, base_name, cmap):
    """
    Plot AUROC heatmap for a given dataset.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing AUROC values indexed by dataset names.
    cell_read_data : array-like
        Shape (n, 2). Each row = (cells_pct, reads_pct).
    dataset_name : str
        Name of the dataset in df.index, e.g. "GENIE3".
    base_name : str
        Prefix for saving the output file.
    cmap : Colormap
        Colormap for the heatmap.
    """

    # Extract values for the given dataset
    values = df.loc[dataset_name].values

    # Build DataFrame
    plot_df = pd.DataFrame(cell_read_data, columns=["cells_pct", "reads_pct"])
    plot_df["cells_pct"] = plot_df["cells_pct"].astype(int)
    plot_df["reads_pct"] = plot_df["reads_pct"].astype(int)
    plot_df["AUROC"] = values

    # Pivot to heatmap format
    heatmap_data = (
        plot_df
        .pivot(index="cells_pct", columns="reads_pct", values="AUROC")
        .sort_index(ascending=False)  # so that highest %cells is at the top
    )

    # Plot
    ax = sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        linewidths=0.5,
        linecolor='lightgray',
        cbar_kws={'label': 'AUROC'},
        square=True,
        annot_kws={'size': 8, 'weight': 'bold', 'color': 'dimgray'}
    )
    cmap = ax.collections[0].cmap
    norm = ax.collections[0].norm

    for text in ax.texts:
        value = float(text.get_text())
        # Map value to RGBA
        rgba = cmap(norm(value))
        # Compute perceived brightness (luma)
        brightness = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
        # Choose black or white depending on brightness
        text.set_color("#2e2d2c" if brightness > 0.5 else "white")

    plt.title(fr"{dataset_name} AUROC across \%Cells and \%Reads for kim23 data", pad=20)
    plt.xlabel(r"\% Reads")
    plt.ylabel(r"\% Cells")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    set_figure_width()
    plt.tight_layout()

    plot_path = f"../plots/{base_name}_{dataset_name.lower()}.pdf"
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()  # close to avoid overlapping if called repeatedly
    print(f"Saved figure: {plot_path}")


for dataset in ["GENIE3", "PIDC"]:
    plot_auroc_heatmap(df, cell_read_data, dataset, base_name, cmap)
