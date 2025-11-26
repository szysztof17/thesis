#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, METHOD_PALETTE_INDICES, get_hue_palette
import seaborn as sns

# === Load style and config ===
load_custom_style(colours = 'default')

method_color_ids = METHOD_PALETTE_INDICES
method_color_ids['CellOracle_noATAC'] = 3
method_color_ids['CellOracle_withATAC'] = 4

palette = get_hue_palette()
method_list_sorted = sorted(method_color_ids, key=method_color_ids.get)

custom_palette_list = [
    palette[method_color_ids[method] ] 
    for method in method_list_sorted
]

# === File paths ===
base_name = "kim23_buenrostro_combined_auroc"
data_path = f"../data/{base_name}.csv"
plot_path = f"../plots/{base_name}.pdf"

# === Load plot-ready data ===
df = pd.read_csv(data_path, index_col=[0, 1])
gene_groups = df.groupby(level=0, sort=False)
n_groups = len(gene_groups)


ax = sns.barplot(df, 
    x="dataset",
    y="AUROC",
    hue="Algorithm",
    hue_order=method_list_sorted,
    palette=custom_palette_list
    )

# --- Add Title ---
ax.set_title('AUROC on Kim23 and Buenrostro18')

ax.set_ylim(0, 1) 

ax.axhline(0.5, color='gray', linestyle='--', linewidth=0.5, label='Random Baseline')

for container in ax.containers:
    ax.bar_label(container, fmt='%.3f', padding=3, fontsize=8)

ax.set_xlabel('Dataset')
ax.set_ylabel('AUROC ')


# --- NEW: RENAME LEGEND ENTRIES ---
handles, labels = ax.get_legend_handles_labels()

rename_map = {
    'CellOracle_noATAC': 'CellOracle (RNA)',
    'CellOracle_withATAC': 'CellOracle (RNA+ATAC)'
}

new_labels = [rename_map.get(label, label) for label in labels]

ax.legend(handles, new_labels, 
          title='Method',
          bbox_to_anchor=(1.02, 0.5), 
          loc='center left',
          borderaxespad=0.,
          fontsize = 'small'
          )


plt.tight_layout()
set_figure_width()
plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
#plt.show()



