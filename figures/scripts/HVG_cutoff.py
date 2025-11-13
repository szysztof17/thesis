#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette 

# === Load style and config ===
load_custom_style()

# === File paths ===
base_name = "HVG_cutoff"
auroc_df_sorted = pd.read_csv('../data/HVG_cutoff.csv')
plot_path = f"../plots/{base_name}.pdf"

#auroc_df_sorted = pd.read_csv('../data/HVG_cutoff.csv')
sns.lineplot(x='NumberOfGenes', y='AUROC', data=auroc_df_sorted, marker='o', hue='Network') 
plt.title('AUROC for Different Number of Genes')
plt.xlabel('NumberOfGenes')
plt.ylabel('AUROC')
plt.xticks(auroc_df_sorted['NumberOfGenes']) # Opcjonalnie

#plt.ylim(0.5, 1.0)

set_figure_width(aspect_ratio=2.25)
plt.tight_layout()

plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
