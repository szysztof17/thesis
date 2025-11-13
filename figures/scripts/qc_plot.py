
import sys
from pathlib import Path

# Add thesis root directory to sys.path
project_root = Path(__file__).resolve().parents[2]  # assuming: thesis/src/jupyter/qc_plot.py
sys.path.insert(0, str(project_root))

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter

from plot_utils import figures, figure_width, load_custom_style, set_figure_width, get_hue_palette, thousands_formatter 

# === Load style and config ===
load_custom_style()


# HM-1
dataset_name = 'kim23-hm1'
prefix = 'GSM4832483_HM-1_'

# --- 1. Load Data ---
data_dir = '/home/ajank/sc_GRN_reconstruction/data/Kim2023'
adata_rna = sc.read_10x_mtx(data_dir, prefix=prefix)

print(f"Initial AnnData shape: {adata_rna.shape}") # (cells, genes)
min_genes = 1259
min_counts = 4000 if dataset_name == 'kim23-dm1' else 100
min_cells = 10
max_mt = 25
## --- 2. Preprocess Data ---

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata_rna.var["mt"] = adata_rna.var_names.str.startswith("MT-")
# ribosomal genes
adata_rna.var["ribo"] = adata_rna.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata_rna.var["hb"] = adata_rna.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata_rna, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

fig = sc.pl.violin(
    adata_rna,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    show=False 
)
# Add horizontal lines to each subplot
axes = fig.axes.flat
lines = [min_genes, min_counts, max_mt] 

# Add line on the first subplot
axes[0].axhline(min_genes, color='red', linestyle='--', linewidth=2, alpha = 0.5)

# Skip middle plot (axes[1]) — no line here

axes[2].axhline(max_mt, color='red', linestyle='--', linewidth=2, alpha = 0.5)

plt.tight_layout()
plt.show()
## --- 2. Preprocess Data ---

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata_rna.var["mt"] = adata_rna.var_names.str.startswith("MT-")
# ribosomal genes
adata_rna.var["ribo"] = adata_rna.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata_rna.var["hb"] = adata_rna.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata_rna, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

fig = sc.pl.violin(
    adata_rna,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    size=0.5,
    alpha=0.5,

    show=False 
)

# Add horizontal lines to each subplot
axes = fig.axes.flat


lines = [min_genes, min_counts, max_mt] 


# === MARK REJECTED REGIONS ===
# First plot: genes < min_genes
axes[0].axhspan(0, min_genes, color='red', alpha=0.2, clip_on=True)

# Third plot: pct_counts_mt > max_mt
axes[2].axhspan(max_mt, axes[2].get_ylim()[1], color='red', alpha=0.2, clip_on=True)

# Horizontal threshold lines
axes[0].axhline(min_genes, color='red', linestyle='--', linewidth=2, alpha = 0.5)
axes[2].axhline(max_mt, color='red', linestyle='--', linewidth=2, alpha = 0.5)



for ax in axes:
    for artist in ax.get_children():
        if isinstance(artist, (mpl.collections.PathCollection, mpl.collections.PolyCollection)):
            artist.set_rasterized(True)
for ax in axes[1:]:  
    ax.set_ylabel("")

# for ax in axes:
#     labels = ax.get_yticks()
#     ax.set_yticklabels([int(lab / 1000) for lab in labels])
#     ax.set_ylabel("value (× 1000)")


for ax in axes:
    ax.yaxis.set_major_formatter(FuncFormatter(thousands_formatter))

figures = Path('/home/kl467102/thesis/figures')
plot_path = figures / '../plots/QC_plot.pdf'
plt.tight_layout()
plt.subplots_adjust(wspace=0.2)  # spacing between subplots

set_figure_width(aspect_ratio=3)
plt.savefig(plot_path, bbox_inches='tight')
print(f"Save a figure to: {plot_path}")
plt.show()
