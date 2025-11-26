#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D  
import seaborn as sns
from pathlib import Path
import os
import re
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette, METHOD_PALETTE_INDICES, get_method_color

import numpy as np
import re

# === Load style and config ===
load_custom_style()

raw_palette = get_hue_palette(style = 'two') 
METHOD_PALETTE = {
    "greedy": raw_palette[0], 
    "random": raw_palette[1]
}

base_name = 'kim23_greedy_random'
# 1. Load Data
df = pd.read_csv(f'../data/{base_name}.csv', index_col=0)
plot_path = f"../plots/{base_name}.pdf"


ax = sns.lineplot(
    data=df,
    x="umi_pct",
    y="AUROC",
    hue="method",
    palette=METHOD_PALETTE,
    marker="o",
    markersize=8,
    linewidth=2.5
)

# for _, row in df.iterrows():
#     method = row["method"]
#     n_cells = int(row["n_cells"])
    
#     if method == "greedy":
#         xytext = (0, 10)   # Shift up
#         va = "bottom"
#     else:
#         xytext = (0, -15)  # Shift down
#         va = "top"         

#     ax.annotate(
#         f"{n_cells}", 
#         (row["umi_pct"], row["AUROC"]),
#         textcoords="offset points",
#         xytext=xytext,
#         ha='center',
#         va=va,
# #        fontsize=9,
#         fontweight='bold',
#         color=METHOD_PALETTE[method], # Match text color to line color
#     )


unique_pcts = df['umi_pct'].unique()

for x in unique_pcts:
    # Get the data for this specific X value
    subset = df[df['umi_pct'] == x]
    
    # We need to determine which method is physically higher at this point
    # separate the rows
    greedy_row = subset[subset['method'] == 'greedy']
    random_row = subset[subset['method'] == 'random']
    
    # Process the rows (handle case where one might be missing, though unlikely)
    current_rows = [greedy_row, random_row]
    
    for row_data in current_rows:
        if row_data.empty: continue
        
        # Extract series to object
        row = row_data.iloc[0] 
        method = row["method"]
        n_cells = int(row["n_cells"])
        y_val = row["AUROC"]
        
        # --- DYNAMIC POSITIONING LOGIC ---
        # If both exist, compare their Y values to decide direction
        if not greedy_row.empty and not random_row.empty:
            greedy_y = greedy_row.iloc[0]['AUROC']
            random_y = random_row.iloc[0]['AUROC']
            
            # If this specific row is the HIGHER one, push text UP
            # If this specific row is the LOWER one, push text DOWN
            if y_val >= max(greedy_y, random_y):
                xytext = (0, 10)
                va = "bottom"
            else:
                xytext = (0, -8)
                va = "top"
        else:
            # Fallback if only one point exists at this X: default to method logic
            if method == "greedy":
                xytext = (0, 10); va = "bottom"
            else:
                xytext = (0, -15); va = "top"
        # ---------------------------------

        ax.annotate(
            f"{n_cells}", 
            (row["umi_pct"], row["AUROC"]),
            textcoords="offset points",
            xytext=xytext,
            ha='center',
            va=va,
            fontweight='bold',
            color=METHOD_PALETTE[method],
            # Adding a white outline makes it readable even if it touches a line
            path_effects=[path_effects.withStroke(linewidth=2, foreground="white")]
        )



handles, labels = ax.get_legend_handles_labels()
info_handle = Line2D([0], [0], color='none', label='Numbers indicate cell count (n)')
handles.append(info_handle)
ax.legend(handles=handles, title="Cells selection", loc="best", fontsize='medium')


plt.title("PIDC AUROC: Greedy vs random sampling evaluated on Kim1000 net")
plt.xlabel(r"\% UMIs (Read depth)")
plt.ylabel("AUROC")
plt.ylim((0.50, 0.85))
plt.grid(True, linestyle='--', alpha=0.5)



set_figure_width(aspect_ratio = 2)
plt.tight_layout()   # reserve right margin for legend
plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
