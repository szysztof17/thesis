#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette 
import numpy as np



# === Load style and config ===
load_custom_style(colours = 'one')

# === File paths ===
base_name = "cross_dataset_auroc"
df = pd.read_csv(f'../data/{base_name}.csv')
plot_path = f"../plots/{base_name}.pdf"

# Plot configuration

datasets = ['kim23-hm1', 'hhep']
tf_sizes = ['500TFs', '1000TFs']
methods = df.index.tolist()

# Map evaluation cases for each dataset and TF size
def get_columns(eval_dataset, tf_size):
    if eval_dataset == 'kim23-hm1':
        same = f"kim23-hm1_STRING12_{tf_size}"
        cross = f"kim23-hm1_hhep{tf_size.replace('TFs', 'Network_' + tf_size)}"
    else:
        same = f"hhep_STRING12_{tf_size}"
        cross = f"hhep_kim{tf_size.replace('TFs', 'Network_' + tf_size)}"
    return same, cross

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for ax, tf_size in zip(axes, tf_sizes):
    x = np.arange(len(datasets))
    width = 0.35

    same_medians, same_errors = [], []
    cross_medians, cross_errors = [], []

    for ds in datasets:
        same_col, cross_col = get_columns(ds, tf_size)

        same_vals = df[same_col]
        cross_vals = df[cross_col]

        same_medians.append(same_vals.median())
        same_errors.append([[same_vals.median() - same_vals.min()], [same_vals.max() - same_vals.median()]])
        cross_medians.append(cross_vals.median())
        cross_errors.append([[cross_vals.median() - cross_vals.min()], [cross_vals.max() - cross_vals.median()]])

    same_errors = np.array(same_errors).squeeze().T
    cross_errors = np.array(cross_errors).squeeze().T

    ax.bar(x - width/2, same_medians, width, yerr=same_errors, capsize=5, label='Dataset specific network')
    ax.bar(x + width/2, cross_medians, width, yerr=cross_errors, capsize=5, label='Cross network')

    ax.set_xticks(x)
    ax.set_xticklabels(datasets)
    ax.set_title(tf_size.replace("TFs", " + TFs"))
    ax.set_ylabel("AUROC")
    ax.set_ylim(0.5, 1)

axes[0].legend()
plt.suptitle("AUROC on in-Dataset vs cross-dataset networks")
plt.tight_layout()
plt.show()







set_figure_width(aspect_ratio=2.25)
plt.tight_layout()

plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
