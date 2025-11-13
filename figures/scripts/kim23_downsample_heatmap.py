#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path
from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette 
import numpy as np


# === Load style and config ===
load_custom_style()

palette = get_hue_palette(6)

# === File paths ===
base_name = "kim23_downsample_heatmap"
df_final = pd.read_csv(f'../data/kim23_budget_constrained.csv')



set_figure_width()
plt.tight_layout()
plot_path = f"../plots/{base_name}_sub.pdf"
plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")


