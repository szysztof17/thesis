#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette 
import numpy as np
from matplotlib.ticker import PercentFormatter

# === Load style and config ===
load_custom_style()

palette = get_hue_palette(6)

# === File paths ===
base_name = "kim23_budget_constrained"
df_final = pd.read_csv(f'../data/kim23_budget_constrained.csv')

mean_auroc_df = df_final.groupby(['budget_name', 'strategy', 'n_cells'])['auroc'].mean().reset_index()

# Set the budget order for plotting
budget_order = ['Nano10', 'Low25', 'Medium50', 'High75']
mean_auroc_df['budget_name'] = pd.Categorical(mean_auroc_df['budget_name'], categories=budget_order, ordered=True)

method_order = ['PIDC', 'GENIE3', 'GRNBOOST2']
df_final['method'] = pd.Categorical(df_final['method'], categories=method_order, ordered=True)

# Sort the data for clean line plotting
mean_auroc_df = mean_auroc_df.sort_values(['budget_name', 'strategy'])

df = df_final.copy()
mean_umis = df.groupby('budget_name')['total_umis_in_data'].mean().to_dict()
budget_order = ['Nano10', 'Low25', 'Medium50', 'High75']
df['budget_name'] = pd.Categorical(df['budget_name'], categories=budget_order, ordered=True)
df = df.sort_values(['budget_name', 'strategy', 'method'])


sns.lineplot(
    data=mean_auroc_df,
    x='strategy',
    y='auroc',
    hue='budget_name',
    style='budget_name', # Use different line styles for each budget
    markers=True,        # Mark each data point
    palette=palette,   # A color-blind friendly palette
    linewidth=2
)
plt.title('Mean performance vs. sequencing strategy')
plt.xlabel('Strategy (0 = Depth-first, 100 = Cells-first)')
plt.ylabel('Mean AUROC')
plt.legend(title='Budget')
set_figure_width(aspect_ratio = 2)
plt.tight_layout()
plot_path = f"../plots/{base_name}_mean.pdf"
plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")



g = sns.relplot(
    data=df,
    x='strategy',
    y='auroc',
    hue='method',
    col='budget_name',  # This creates the subplots for each budget
    kind='line',        # Specify the plot type
    col_wrap=2,         # Wrap the subplots into a 2x2 grid
    palette=palette,    # Color palette for the lines
    marker='o',         # Style for the markers
    height=4,           # Height of each subplot
    aspect=1.2,         # Aspect ratio of each subplot
    linewidth=1,
    markersize = 4
)
for ax in g.axes.flatten():
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
g.fig.suptitle('Method performance vs. sequencing strategy across budgets')#, y=1.03)
g.set_axis_labels("", "AUROC") # Clear individual labels first
g.fig.text(0.5, -0.01, r'Strategy (0\% = Depth-first, 100\% = Cells-first)', ha='center', va='bottom')

g._legend.set_bbox_to_anchor((1, 0.2))

for budget, ax in g.axes_dict.items():
    avg_umis = mean_umis.get(budget, 0)
    title = rf"Budget: {budget},  $\sim${avg_umis/1e6:.1f}M UMIs"
    ax.set_title(title)


set_figure_width()
plt.tight_layout()
plot_path = f"../plots/{base_name}_sub.pdf"
plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")


