import os
import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
from tqdm import tqdm

# --- Helper Functions (from step3_downsampling.ipynb) ---

def sample_cells(adata, *, n_cells=None, frac=None, seed=420):
    """
    Returns a copy of `adata` containing a random subset of cells.
   
    """
    if (n_cells is None) == (frac is None):
        raise ValueError("specify *either* n_cells or frac, not both")

    rng = np.random.default_rng(seed)
    total = adata.n_obs

    if frac is not None:
        if not (0 < frac <= 1):
            raise ValueError("`frac` must be in (0, 1]")
        n_cells = int(round(total * frac))

    if n_cells > total:
        print(f"Warning: requested {n_cells} cells but adata only has {total}. Using all {total} cells.")
        n_cells = total
    
    keep = rng.choice(total, size=n_cells, replace=False)
    return adata[keep].copy()

    
def sample_reads(adata, *, depth_factor=None, target_umi=None, seed=420):
    """
    Returns a copy of `adata` with down-sampled counts.
   
    """
    if (depth_factor is None) == (target_umi is None):
        raise ValueError("specify exactly one of depth_factor or target_umi")

    rng = np.random.default_rng(seed)
    X = adata.layers['counts'].copy()

    if sp.issparse(X): # handle sparse matrices efficiently
        X = X.tocoo()
        if depth_factor is not None:
            p = depth_factor
            X.data = rng.binomial(X.data.astype(np.int64), p)
        else:
            totals = np.asarray(adata.X.sum(axis=1)).ravel() # per-cell UMIs
            p_cell = np.minimum(1., target_umi / totals) # per-cell probs
            X.data = rng.binomial(X.data.astype(np.int64), p_cell[X.row])
        X = X.tocsr()
    else: # dense matrix
        if depth_factor is not None:
            X = rng.binomial(X.astype(np.int64), depth_factor)
        else:
            totals = X.sum(1, keepdims=True)
            p_cell = np.minimum(1., target_umi / totals)
            X = rng.binomial(X.astype(np.int64), p_cell)

    out = adata.copy()
    out.X = X
    out.obs["n_counts_sim"] = np.asarray(X.sum(axis=1)).ravel() # keep track
    return out

def build_hvg_subset(ad, n_top_hvg, all_known_tfs, include_hv_tfs):
    """
    Identify & combine HVGs (from step3_downsampling.ipynb).
   
    """
    adata = ad.copy()
    
    # Normalize and logarithmize
    adata.raw = adata.copy()
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Default HVGs
    sc.pp.highly_variable_genes(adata, n_top_genes=None, subset=False)
    default_mask = adata.var['highly_variable'].copy()
    n_default = default_mask.sum()
    print(f"Default HVGs: {n_default}")
    
    # Cutoff HVGs
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_hvg, subset=False)
    cutoff_mask = adata.var['highly_variable'].copy()
    n_cutoff = cutoff_mask.sum()
    print(f"Top {n_top_hvg} HVGs: {n_cutoff}")
    
    # Include TFs
    tf_mask = adata.var_names.isin(all_known_tfs)
    tf_in_default = default_mask & tf_mask
    n_tf_default = tf_in_default.sum()
    print(f"TFs in default HVGs: {n_tf_default}")
    
    combined = cutoff_mask.copy()
    if include_hv_tfs:
        combined |= tf_in_default
        added = combined.sum() - n_cutoff
        print(f"Including additional {added} TFs → total combined HVGs: {combined.sum()}")
        
    ad_hvg = adata[:, combined].copy()
    # ...
    if sp.issparse(ad_hvg.X):
        X = ad_hvg.X.toarray()
    else:
        X = np.array(ad_hvg.X)

    expr_df = pd.DataFrame(
        X,
        index=ad_hvg.obs_names,
        columns=ad_hvg.var_names
    ).T
    print(f"HVG subset shape: cells={ad_hvg.n_obs}, genes={ad_hvg.n_vars}\n")
    return ad_hvg, expr_df


def filter_and_save(ad_sub, expr_df, net_df, cfg, key):
    """
    Filter & save network + expression (from step3_downsampling.ipynb).
   
    """
    print("--- STEP 3: Filtering network & saving outputs ---")
    cells_frac, reads_frac = key
    
    # Dataset part with hyphens: DATASET-CELLS-READS
    dataset_part = f"{cfg['dataset_name']}-{cells_frac}-{reads_frac}"
    
    # Suffix: _NETWORK_NUMGENES(-TFs)
    tfs_suffix = "TFs" if cfg['include_hv_tfs'] else ""
    suffix = f"_{cfg['network_name']}_{cfg['n_top_hvg']}{tfs_suffix}"
    
    prefix = cfg['output_dir'] / f"{dataset_part}{suffix}"
    prefix.mkdir(parents=True, exist_ok=True)
    

    # Filter network
    filtered_net_df = net_df[
        net_df['Gene1'].isin(ad_sub.var.index) &
        net_df['Gene2'].isin(ad_sub.var.index)
    ].copy()
    filtered_net_df = filtered_net_df[filtered_net_df['Gene1'] != filtered_net_df['Gene2']]
    filtered_net_df.drop_duplicates(keep='first', inplace=True)
    
    # Save files
    expr_file = prefix / 'ExpressionData.csv'
    expr_df.to_csv(expr_file, index=True)
    filtered_network_filename = prefix / 'refNetwork.csv'
    filtered_net_df.to_csv(filtered_network_filename, index=False)
    adata_storade = Path('adatas_to_be_deleted_later')
    os.makedirs(adata_storade, exist_ok=True)
    h5ad_file =  adata_storade / f'{prefix}.h5ad'
    ad_sub.write(h5ad_file)


    print(f"{cells_frac},{reads_frac} → cells:{ad_sub.n_obs}, genes:{ad_sub.n_vars}")


    return prefix.name





def main(args):
    # --- Config ---
    cfg = {
        'tf_file': args.tf_file,
        'network_file': args.network_file,
        'network_name': args.network_name,
        'output_dir': Path(args.beeline_dir) / 'inputs' / args.dataset_name,
        'dataset_name': args.dataset_name,
        'cell_fracs': np.linspace(0.1, 1.0, 10), #
        'read_fracs': np.linspace(0.1, 1.0, 10), #
        'seed': 420,
        'n_top_hvg': args.n_top_hvg,
        'include_hv_tfs': args.include_hv_tfs,
    }
    cfg['output_dir'].mkdir(parents=True, exist_ok=True)

    # --- Load global files ---
    print(f"Reading TF list: {cfg['tf_file']}")
    all_known_tfs = pd.read_csv(cfg['tf_file'], header=None).iloc[:, 0].tolist()
    print(f"Reading network: {cfg['network_file']}")
    net_df = pd.read_csv(cfg['network_file'])
    
    # --- Load main dataset ---
    print(f"Loading data from: {args.input_data}")
    # This assumes the new data is in the same .txt format as buenrostro18
    # If loading .h5ad, use sc.read() instead.
    adata_rna =  sc.read(args.input_data)
    print(f"Initial AnnData shape: {adata_rna.shape}")

    # --- Run Pipeline ---
    datasets = []
    print("\n--- STEP 1: Subsampling reads over grid ---")
    for f_cells in tqdm(cfg['cell_fracs'], desc="Sampling cell fractions"):
        ad_sub = sample_cells(adata_rna, frac=f_cells, seed=cfg['seed'])
        
        for f_reads in tqdm(cfg['read_fracs'], desc=f"Reads for {int(f_cells*100)}% cells", leave=False):
            ad_sim = sample_reads(ad_sub, depth_factor=f_reads, seed=cfg['seed'])
            
            # STEP 2: HVG selection
            ad_hvg, expr_df = build_hvg_subset(
                ad_sim, cfg["n_top_hvg"], all_known_tfs, cfg["include_hv_tfs"]
            )
            
            # STEP 3: Filter and save
            key = (f"{int(f_cells*100)}pct-cells", f"{int(f_reads*100)}pct-reads")
            dataset_name = filter_and_save(ad_hvg, expr_df, net_df, cfg, key)
            datasets.append(dataset_name)

    print("\n✅ Pipeline finished.")
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modular BEELINE preprocessing and downsampling script.")
    
    # --- Required Paths ---
    parser.add_argument('--input_data', type=str, required=True, 
                        help="Path to the input .h5ad")
    parser.add_argument('--dataset_name', type=str, required=True, 
                        help="Name for the new dataset (e.g., 'new_data_2025')")

    # --- Paths with defaults (change if needed) ---
    parser.add_argument('--tf_file', type=str, 
                        default='/home/kl467102/Beeline-238/inputs/hHep/human-tfs.csv', 
                        help="Path to human TFs list")
    parser.add_argument('--network_file', type=str, 
                        default='/home/kl467102/string_dir/9606_protein_links_gene_names_combined_score_700.csv', 
                        help="Path to reference network")
    parser.add_argument('--beeline_dir', type=str, 
                        default='/home/kl467102/thesis/BEELINE/', 
                        help="BEELINE directory path")

    # --- Parameters ---
    parser.add_argument('--n_top_hvg', type=int, default=1000, 
                        help="Number of HVGs")
    parser.add_argument('--network_name', type=str, default='STRING12', 
                        help="Network name suffix")
    
    # --- Flags ---
    parser.add_argument('--no_include_hv_tfs', dest='include_hv_tfs', action='store_false',
                        help="Flag to NOT include TFs from default HVG list")
    parser.set_defaults(include_hv_tfs=True)

    
    args = parser.parse_args()
    main(args)