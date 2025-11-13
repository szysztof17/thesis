#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width

# === Load style and config ===
load_custom_style()

# === File paths ===
base_name = "hhep_methods_comparison"
data_path = f"../data/{base_name}1.csv"
plot_path = f"../plots/{base_name}.pdf"

# === Load plot-ready data ===
df = pd.read_csv(data_path, index_col=[0, 1])
df.sort_index(level=0, key=lambda x: [int(s.replace("TFs", "")) for s in x], inplace=True)
gene_groups = df.groupby(level=0, sort=False)
n_groups = len(gene_groups)

'''
fig, axes = plt.subplots(1, n_groups, sharey=True)

if n_groups == 1:
    axes = [axes]

for ax, (gene_count, group) in zip(axes, gene_groups):
    x = range(len(group))
    width = 0.25

    ax.bar([i - width for i in x], group['STRING12'], width=width, label='STRING12')
    ax.bar(x, group['universalNet'], width=width, label='universalNet')
    ax.bar([i + width for i in x], group['minimalNet'], width=width, label='minimalNet')

    ax.set_title(f'Number of genes: {gene_count.replace('TFs','')} + TFs')
    ax.set_xticks(list(x))
    ax.set_xticklabels(group.index.get_level_values('Algorithm'), rotation=45, ha='right')
    ax.set_ylim(0, 1)
    ax.set_ylabel('AUROC')
    ax.legend()
'''
fig, axes = plt.subplots(1, n_groups, sharey=True)

# Ensure axes is iterable
if n_groups == 1:
    axes = [axes]
i=0
for ax, (gene_count, group) in zip(axes, gene_groups):
    i+=1
    x = range(len(group))
    width = 0.25

    ax.bar([i - width for i in x], group['hhep1000Net'], width=width, label='1000+TFs Net')
    ax.bar(x, group['hhep500Net'], width=width, label='500+TFs Net')
    ax.bar([i + width for i in x], group['minimalNet'], width=width, label='minimalNet')

    ax.set_title(f'Number of genes: {gene_count.replace('TFs','')} + TFs ({int(gene_count.replace('TFs','')) + 448})')
    ax.set_xticks(list(x))
    ax.set_xticklabels(group.index.get_level_values('Algorithm'), rotation=45, ha='right')
    ax.set_ylim(0, 1)
    ax.set_ylabel('AUROC')
    if i == n_groups:
        ax.legend(title="Reference Net.")

plt.tight_layout()
set_figure_width(aspect_ratio=2)
plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
#plt.show()
