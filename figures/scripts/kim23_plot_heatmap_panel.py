#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
from pathlib import Path
import re


label_size = 8
ticks_size = 7
title_size = 9
data_threshold = 0.76
try:
    from plot_utils import load_custom_style, get_method_color, make_sequential_cmap, set_figure_width
except ImportError:
    print("ERROR: Could not import from plot_utils.py.")
    print("Please ensure plot_utils.py is in the same directory and configured correctly.")
    exit()

# === 1. Helper Plotting Function (Modified from your original) ===

def plot_heatmap_panel(ax, df, cell_read_data, dataset_name, cmap):
    """
    Plots a *single* heatmap panel onto a given matplotlib Axes (ax).
    Returns the mappable object needed to create a shared colorbar.
    """
    # --- Data Preparation ---
    # Assumes 'df' is indexed by method name ('dataset_name')
    values = df.loc[dataset_name].values
    plot_df = pd.DataFrame(cell_read_data, columns=["cells_pct", "reads_pct"])
    plot_df["cells_pct"] = plot_df["cells_pct"].astype(int)
    plot_df["reads_pct"] = plot_df["reads_pct"].astype(int)
    plot_df["AUROC"] = values
    
    # Pivot for heatmap
    heatmap_data = plot_df.pivot(index="cells_pct", columns="reads_pct", values="AUROC").sort_index(ascending=False)

    # --- Plot Heatmap ---
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        linewidths=0.5,
        cbar=False,  # --- IMPORTANT: Colorbar is handled externally ---
        ax=ax,
        annot_kws={'size': 5} # Smaller annotations for a dense plot
    )
    
    # --- Text Color Logic (Using the data threshold) ---
    # This ensures text color flips at the same AUROC value for all plots
    
    for text in ax.texts:
        value = float(text.get_text())
        text.set_color("#2e2d2c" if value < data_threshold else "white")
        # Format without leading zero, 2 decimal places
        text.set_text(f"{value:.2f}".lstrip("0"))
        text.set_color("#2e2d2c" if value < data_threshold else "white")

    ax.set_title(dataset_name, fontsize = title_size)
    ax.set_xlabel(r"\% Reads", fontsize = label_size)
    ax.set_ylabel(r"\% Cells", fontsize = label_size)
    ax.tick_params(axis='both', labelsize = ticks_size)

    
    # Return the mappable object (the heatmap artist)
    # This is used to draw the shared colorbar later
    return ax.collections[0]

# === 2. Main Script ===

def main():
    # --- Load Style ---
    load_custom_style()
    
    # --- Configuration ---
    base_name = "../data/kim23_heatmap_panel.csv"
    methods = ["PIDC", "GRNBOOST2", "GENIE3"]
    
    # === 3. Load Your Data ===
    df = pd.read_csv('../data/kim23_downsample_heatmap.csv', index_col = 0, header = 0)
    df.index.name = "dataset_name"
    dataset_names = df.columns.tolist()
    pattern = re.compile(r"(\d+)pct-cells-(\d+)pct-reads")
    cell_read_data = [pattern.search(name).groups() for name in dataset_names]



    
    # === 4. Figure Setup ===
    
    fig, axes = plt.subplots(1, 3)
    
    mappable = None 
    
    # === 5. Plotting Loop ===
    print("Generating plots...")
    for i, method in enumerate(methods):
        ax = axes[i]
        # Get consistent color and cmap from plot_utils
        try:
            color = get_method_color(method)
            cmap = make_sequential_cmap(color)
        except Exception as e:
            print(f"Warning: Could not get method color for '{method}'. Defaulting. Error: {e}")
            cmap = "viridis" # Fallback
            
        # Plot the panel
        mappable = plot_heatmap_panel(
            ax, df, cell_read_data, method, cmap        
        )
        
        # Clean up Y-axis labels for 2nd and 3rd plots to reduce clutter
        if i > 0:
            ax.set_ylabel("")
            

    fig.suptitle("AUROC Consistency Across Methods (Downsampling)")

    set_figure_width(aspect_ratio=2.25)
    plt.tight_layout()

    plot_path = f"../plots/{base_name}.pdf"
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Saved comparison figure to: {plot_path}")

if __name__ == "__main__":
    main()