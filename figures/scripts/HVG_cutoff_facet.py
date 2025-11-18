#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_utils import load_custom_style, set_figure_width

# === Load style ===
load_custom_style(colours = 'one')

# === File paths ===
csv_path = "../data/HVG_cutoff_kim_hhep.csv"
plot_path = "../plots/HVG_cutoff_facet.pdf"

# === Load Data ===
df = pd.read_csv(csv_path)
df['Network'] = df['Network'].replace({'commonnetwork': 'hhep500Net', 'STRING12': 'Dataset-specific'})
df_hhep = df[df['Dataset']=='hhep'].sort_values(by = "Network")
df_kim = df[df['Dataset']=='kim23-hm1'].sort_values(by = "Network")
# === Plotting ===
fig, axes = plt.subplots(1, 2)  

sns.lineplot(x='NumberOfGenes', y='AUROC', data=df_hhep, marker='o', hue='Network', ax = axes[0]) 
axes[0].set_title("hhep")
axes[0].set_xlabel("Number of genes")
axes[0].set_ylabel("AUROC")
axes[0].legend(title=r"Ref. network", fontsize='small')



sns.lineplot(x='NumberOfGenes', y='AUROC', data=df_kim, marker='o', hue='Network', ax = axes[1]) 
axes[1].set_title("kim23")
axes[1].set_xlabel("Number of genes")
axes[1].set_ylabel("AUROC")
axes[1].legend(title=r"Ref. network", fontsize='small')


axes[0].set_ylim(0.5, 1.0)
axes[1].set_ylim(0.5, 1.0)

axes[1].set_ylabel("")
axes[1].tick_params(labelleft=False)

set_figure_width(aspect_ratio=2.25)
plt.tight_layout()

plt.savefig(plot_path, bbox_inches='tight')
print(f"Saved figure to: {plot_path}")
