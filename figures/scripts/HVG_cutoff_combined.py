#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_utils import load_custom_style, set_figure_width

# === Load style ===
load_custom_style()

# === File paths ===
csv_path = "../data/HVG_cutoff_kim_hhep.csv"
plot_path = "../plots/HVG_cutoff_combined.pdf"

# === Load Data ===
df = pd.read_csv(csv_path)

# Create combined label for color legend
df["Label"] = df["Dataset"] + ":" + df["Network"]

# === Plotting ===
set_figure_width(aspect_ratio=2.25)
sns.lineplot(
    data=df,
    x="NumberOfGenes",
    y="AUROC",
    hue="Label",
    marker="o"
)

plt.title("AUROC vs Number of Genes")
plt.xlabel("Number of Genes")
plt.ylabel("AUROC")
plt.xticks(sorted(df["NumberOfGenes"].unique()))
plt.tight_layout()

plt.savefig(plot_path, bbox_inches='tight')
print(f"Saved figure to: {plot_path}")
