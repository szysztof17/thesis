#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette, make_sequential_cmap
import numpy as np
import re

# === Load style and config ===
load_custom_style()

palette = get_hue_palette(6)

cmap = make_sequential_cmap(sns.color_palette()[1])
#cmap = make_cmap_between(palette[0], palette[2])


# === File paths ===
base_name = "kim23_meta_grid_user_times"
df = pd.read_csv(f'../data/{base_name}.csv', index_col=0)
df = df.loc[['GRNBOOST2']]
sns.scatterplot(
    data=df,
    x="user_time_seconds",
    y="auroc",
    hue="method",
    style="metacells",
    palette = palette[2:],
    s=20  # marker size
)
plt.xlabel("Running time (seconds)")
plt.ylabel("AUROC")
plt.title("Trade-off between running time and AUROC")
plt.grid(True)


set_figure_width(n_per_page=2, aspect_ratio = 1)
plt.tight_layout()

plot_path = f"../plots/{base_name}.pdf"
plt.savefig(plot_path, bbox_inches='tight')
print(f"Saved figure: {plot_path}")
