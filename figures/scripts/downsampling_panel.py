#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
from pathlib import Path
import re
import matplotlib.lines as mlines

label_size = 8
ticks_size = 7
title_size = 12
data_threshold = 0.76
try:
    from plot_utils import load_custom_style, get_method_color, make_sequential_cmap, set_figure_width, get_hue_palette, METHOD_PALETTE_INDICES, get_method_color
except ImportError:
    print("ERROR: Could not import from plot_utils.py.")
    print("Please ensure plot_utils.py is in the same directory and configured correctly.")
    exit()


FILES = {
    "Kim23": "../data/kim23_downsample_heatmap.csv",          
    "Buenrostro18": "../data/buenrostro18-auroc.csv",        
    #"PBMC10k": "../data/pbmc10k_downsample_netowrk_AUROC.csv"    
}   
plot_path = "../plots/grid_heatmap_comparison.pdf"
METHODS = ["PIDC", "GRNBOOST2", "GENIE3"]
PATTERN = re.compile(r"(\d+)pct-cells-(\d+)pct-reads")

# Visual Settings
CMAPS  = {m: get_method_color(m) for m in METHOD_PALETTE_INDICES}

VMIN, VMAX = 0.5, 0.8  # Fixed range for fair comparison across ALL plots

def parse_and_format(df, dataset_prefix, method_name):
    """
    Extracts downsampling columns for a specific method and formats into a matrix.
    """
    # 1. Filter columns based on pattern 
    cols = [c for c in df.columns if "pct-cells" in c and PATTERN.search(c)]
    
    if not cols:
        print(f"Warning: No downsampling columns found in {dataset_prefix}")
        return None
        
    if method_name not in df.index:
        print(f"Warning: Method {method_name} not found in index of {dataset_prefix}")
        return None
    # 2. Extract Data
    data = []
    for col in cols:
        match = PATTERN.search(col)
        if match:
            cells, reads = map(int, match.groups())
            auroc = df.loc[method_name, col] if method_name in df.index else df[col].mean()
            data.append({"cells": cells, "reads": reads, "AUROC": auroc})
            
    if not data: return None

    # 3. Pivot to Matrix (Heatmap format)

    matrix = pd.DataFrame(data).pivot(index="cells", columns="reads", values="AUROC")
    matrix = matrix.sort_index(ascending=False) 
    return matrix

def main():
    load_custom_style()
    n_rows, n_cols = 2, 3  

    fig, axes = plt.subplots( n_rows, n_cols, sharex=True, sharey=True)


    for row_idx, (dataset_name, filepath) in enumerate(FILES.items()):
        try:
            df = pd.read_csv(filepath, index_col=0)
        except FileNotFoundError:
            print(f"Skipping {dataset_name} (File not found)")
            continue

        for col_idx, method in enumerate(METHODS):
            ax = axes[row_idx, col_idx]
            
            # Prepare Data
            matrix = parse_and_format(df, dataset_name.lower(), method).round(2)
            color = get_method_color(method)
            cmap = make_sequential_cmap(color)
            local_max = matrix.max().max()
            if matrix is not None:
                # Plot Heatmap
                sns.heatmap(
                    matrix, 
                    annot=True,
                    fmt=".2f",
                    cmap=cmap,
                    linewidths=0.5,
                    cbar=False,  
                    ax=ax,
                    annot_kws={'size': 5, } # Smaller annotations for a dense plot
                )
                for text in ax.texts:
                    value = float(text.get_text())
                    # Format without leading zero, 2 decimal places
                    text.set_text(f"{value:.2f}".lstrip("0"))
                    text.set_color("#2e2d2c" if value < local_max - 0.02 else "white")
                    text.set_weight("normal" if value < local_max - 0.02 else "bold")

            # --- Labeling Logic ---
            
            # Column Titles (Methods) - Only on top row
            if row_idx == 0:
                ax.set_title(method, fontsize = label_size)
            if col_idx == 0:
                ax.set_ylabel(f"{dataset_name}\n\n \\% Cells", fontsize=label_size)
            else:
                ax.set_ylabel("")
            # X-axis Labels - Only on bottom row
            if row_idx == n_rows - 1:
                ax.set_xlabel(r"\% Reads", fontsize = label_size)
            else:
                ax.set_xlabel("")
            ax.tick_params(axis='both', labelsize = ticks_size)





    # for row_idx in range(n_rows - 1):  
    #     ax = axes[row_idx, 0]  
    #     pos = ax.get_position()
    #     y = pos.y0    

    #     line = mlines.Line2D(
    #         [0.05, 0.98],   
    #         [y - 0.02, y - 0.02],         
    #         color='gray', linewidth=1, linestyle='--', transform=fig.transFigure, zorder=20
    #     )
    #     fig.add_artist(line)

    fig.suptitle("AUROC Consistency across methods and datasets (downsampling)", fontsize  = title_size, y=0.95)


    set_figure_width(aspect_ratio=3/2.25)
    plt.tight_layout()

    plt.savefig(plot_path, bbox_inches = 'tight')
    print("Plot saved!")

if __name__ == "__main__":
    main()