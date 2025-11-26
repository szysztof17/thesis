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
method_color_ids['pySCENIC_proteins'] = 5
method_color_ids['SCENICPLUS_new_params'] = 6
method_color_ids['LINGER'] = 7


palette = get_hue_palette()
method_list_sorted = sorted(method_color_ids, key=method_color_ids.get)

custom_palette_list = [
    palette[method_color_ids[method] ] 
    for method in method_list_sorted
]

# === File paths ===
base_name = "pbmc10_full_metrics"
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
ax.set_title('AUROC on PBMC10k')

ax.set_ylim(0, 1) 

for container in ax.containers:
    ax.bar_label(container, fmt='%.3f', padding=1, fontsize=6)

#ax.set_ylabel('AUROC ')

ax.set(xticklabels=[])
ax.set(xlabel=None)
ax.set(ylabel=None)
print('ylabel=None')
# --- NEW: RENAME LEGEND ENTRIES ---
handles, labels = ax.get_legend_handles_labels()

rename_map = {
    'CellOracle_noATAC': 'CellOracle',
    'CellOracle_withATAC': 'CellOracle (+ATAC)',
    'SCENICPLUS_new_params': 'SCENIC+',
    'pySCENIC_proteins' : 'SCENIC'
}

new_labels = [rename_map.get(label, label) for label in labels]

ax.legend(handles, new_labels, 
          title='Method',
          borderaxespad=0.,
          fontsize = 'x-small'
          )


plt.tight_layout()
set_figure_width(aspect_ratio=0.9, n_per_page=2)
fig = plt.gcf() if ax is None else ax.figure

print(f"current dims in inch are: {fig.get_size_inches()}")
plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
#plt.show()



