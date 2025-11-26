import sys
import os
import argparse
from pathlib import Path
from joblib import Parallel, delayed
from datetime import datetime
import time

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.cluster import KMeans

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

import SEACells



def make_metacells(adata, n_clusters, method='sum', random_state = 420):
    # Compute PCA if not already done
    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')
    X = adata.obsm['X_pca']
    print(f"Calculated  PCA")

    # Cluster cells into n_clusters using KMeans
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state).fit(X)
    labels = kmeans.labels_
    print(f"Done KMeans")

    metacell_exprs = []
    metacell_names = []

    for cluster_id in range(n_clusters):
        idxs = np.where(labels == cluster_id)[0]
        expr_block = adata.X[idxs]
        if method == 'mean':
            meta_expr = expr_block.mean(axis=0)
        elif method == 'sum':
            meta_expr = expr_block.sum(axis=0)
        else:
            raise ValueError("method must be 'mean' or 'sum'")

        if sp.issparse(meta_expr):
            meta_expr = meta_expr.A1
        elif hasattr(meta_expr, "A1"):
            meta_expr = meta_expr.A1
        else:
            meta_expr = np.asarray(meta_expr).ravel()

        metacell_exprs.append(meta_expr)
        metacell_names.append(f"meta_{cluster_id}")

    X_meta = np.vstack(metacell_exprs)
    meta_adata = sc.AnnData(X=X_meta)
    meta_adata.var_names = adata.var_names
    meta_adata.obs_names = metacell_names

    return meta_adata, labels, n_clusters





def main(args):
    # --- Config ---
    cfg = {
        'tf_file': args.tf_file,
        'network_file': args.network_file,
        'network_name': args.network_name,
        'output_dir': Path(args.beeline_dir) / 'inputs' / args.dataset_name,
        'dataset_name': args.dataset_name,
        'k_list': args.k_values,
        'seed': 420,
        'n_top_hvg': args.n_top_hvg,
        'include_hv_tfs': args.include_hv_tfs,
    }
    cfg['output_dir'].mkdir(parents=True, exist_ok=True)
    output_dir = cfg['output_dir']
    dataset_name = cfg['dataset_name']
    network_name = cfg['network_name']
    n_top_hvg = cfg['n_top_hvg']
    tfs_suffix = 'TFs' if cfg['include_hv_tfs'] else ''
    include_hv_tfs = cfg['include_hv_tfs']
    # --- Load global files ---
    print(f"Reading TF list: {cfg['tf_file']}")
    all_known_tfs = pd.read_csv(cfg['tf_file'], header=None).iloc[:, 0].tolist()
    print(f"Reading network: {cfg['network_file']}")
    net_df = pd.read_csv(cfg['network_file'])
    
    # --- Load main dataset ---
    print(f"Loading data from: {args.input_data}")
    adata_rna = sc.read(args.input_data)
    print(f"Initial AnnData shape: {adata_rna.shape}")


    adata_normalized = adata_rna.copy()
    n_cells = adata_normalized.shape[0]
    
    print("\n--- Running pipeline for non-metacelled (original) data ---")
    param_prefix_orig = output_dir / f"{dataset_name}_{network_name}_{n_top_hvg}{tfs_suffix}"
    os.makedirs(param_prefix_orig, exist_ok=True)
    
    # Copy original AnnData to avoid in-place modification
    adata_orig = adata_normalized.copy()

    sc.pp.normalize_total(adata_orig, target_sum=1e4)
    sc.pp.log1p(adata_orig)
    sc.pp.highly_variable_genes(adata_orig, n_top_genes=None)

    # PROCEED WITH THE BEELINE-LIKE FEAUTRE SELECTION AND NETWORK FILTERING
    default_mask = adata_orig.var['highly_variable'].copy()
    n_default = default_mask.sum()
    print(f"Default HVGs: {n_default}")
    
    # Cutoff HVGs
    sc.pp.highly_variable_genes(adata_orig, n_top_genes=n_top_hvg)
    cutoff_mask = adata_orig.var['highly_variable'].copy()
    n_cutoff = cutoff_mask.sum()
    
    # Include TFs
    tf_mask = adata_orig.var_names.isin(all_known_tfs)
    tf_in_default = default_mask & tf_mask
    n_tf_default = tf_in_default.sum()
    combined = cutoff_mask.copy()
    if include_hv_tfs:
        combined |= tf_in_default
    
    if include_hv_tfs:
        adata_subset_orig = adata_orig[:, combined].copy()
    else:
        adata_subset_orig = adata_orig[:, adata_orig.var['highly_variable']].copy()
    
    # Handle sparse to dense conversion
    if sp.issparse(adata_subset_orig.X):
        expr_df = pd.DataFrame(adata_subset_orig.X.toarray(), index=adata_subset_orig.obs_names, columns=adata_subset_orig.var_names).T
    else:
        expr_df = pd.DataFrame(adata_subset_orig.X, index=adata_subset_orig.obs_names, columns=adata_subset_orig.var_names).T

    # Filter network
    filtered_net_df = net_df[
        net_df['Gene1'].isin(adata_subset_orig.var.index) &
        net_df['Gene2'].isin(adata_subset_orig.var.index)
    ].copy()
    filtered_net_df = filtered_net_df[filtered_net_df['Gene1'] != filtered_net_df['Gene2']]
    filtered_net_df.drop_duplicates(keep='first', inplace=True)
    
    expr_file = param_prefix_orig / 'ExpressionData.csv'
    expr_df.to_csv(expr_file, index=True)
    filtered_network_filename = param_prefix_orig / 'refNetwork.csv'
    filtered_net_df.to_csv(filtered_network_filename, index=False)
    print(f'Saved non-metacelled data to: {param_prefix_orig.name}')
    print("--- Finished non-metacelled pipeline ---")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # --- End of new block ---
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    datasets = []
    k_list = cfg['k_list']
    results_summary = []
    for k in k_list:
        print(f"\n--- Running pipeline with k = {k} ---")
        n = n_cells // k
        param_prefix = output_dir / f"{dataset_name}-KMeans-nmc{n}-avg{k}_{network_name}_{n_top_hvg}{tfs_suffix}"
        os.makedirs(param_prefix, exist_ok=True)
        # Copy original AnnData to avoid in-place modification
        adata = adata_normalized.copy()

        meta_adata, labels, n_clusters = make_metacells(adata, n, method='sum')
        sc.pp.normalize_total(meta_adata, target_sum=1e4)
        sc.pp.log1p(meta_adata)
        sc.pp.highly_variable_genes(meta_adata, n_top_genes=None)

        #PROCEED WITH THE BEELINE-LIKE FEAUTRE SELECTION AND NETWORK FILTERING

        default_mask = adata.var['highly_variable'].copy()
        n_default = default_mask.sum()
        print(f"Default HVGs: {n_default}")
        # Cutoff HVGs
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_hvg)
        cutoff_mask = adata.var['highly_variable'].copy()
        n_cutoff = cutoff_mask.sum()
        # Include TFs
        tf_mask = adata.var_names.isin(all_known_tfs)
        tf_in_default = default_mask & tf_mask
        n_tf_default = tf_in_default.sum()
        combined = cutoff_mask.copy()
        if include_hv_tfs:
            combined |= tf_in_default
        if include_hv_tfs:
            adata_subset = meta_adata[:, combined].copy()
        else:
            adata_subset = meta_adata[:, adata.var['highly_variable']].copy()
        expr_df = pd.DataFrame(adata_subset.X, index=adata_subset.obs_names, columns=adata_subset.var_names).T
        # Filter network

        filtered_net_df = net_df[
            net_df['Gene1'].isin(adata_subset.var.index) &
            net_df['Gene2'].isin(adata_subset.var.index)
        ].copy()
        filtered_net_df = filtered_net_df[filtered_net_df['Gene1'] != filtered_net_df['Gene2']]
        filtered_net_df.drop_duplicates(keep='first', inplace=True)
        expr_file = param_prefix / 'ExpressionData.csv'
        expr_df.to_csv(expr_file, index=True)
        filtered_network_filename = param_prefix / 'refNetwork.csv'
        filtered_net_df.to_csv(filtered_network_filename, index=False)

        datasets.append(param_prefix.name)
        
        results_summary.append({
            'n_cells': adata_subset.n_obs,
            'n_genes': adata_subset.n_vars,
            'k': k,
            'n_clusters': n_clusters,
            'labels': labels,
        })

    # Display summary table
    summary_df = pd.DataFrame(results_summary)
    print(summary_df.head(5))

    
    def run_pipeline(k):
        n_SEACells = n_cells // k
        param_prefix = output_dir / f"{dataset_name}-SEACells-nmc{n_SEACells}-avg{k}_{network_name}_{n_top_hvg}{tfs_suffix}"
        os.makedirs(param_prefix, exist_ok=True)

        start_time = time.time()
        print(f"[{datetime.now().strftime('%H:%M:%S')}] START k={n_SEACells}")

        adata = adata_normalized.copy()

        model = SEACells.core.SEACells(adata,
                        build_kernel_on='X_pca',
                        n_SEACells=n_SEACells,
                        n_waypoint_eigs=10,
                        convergence_epsilon=1e-5)
        
        model.construct_kernel_matrix()
        model.initialize_archetypes()
        model.fit(min_iter=10, max_iter=100)
        for _ in range(5):
            model.step()

        SEACell_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')
        sc.pp.normalize_total(SEACell_ad, target_sum=1e4)
        sc.pp.log1p(SEACell_ad)
        sc.pp.highly_variable_genes(SEACell_ad, n_top_genes=None)
        
        #PROCEED WITH THE BEELINE-LIKE FEAUTRE SELECTION AND NETWORK FILTERING

        adata = SEACell_ad.copy()
        default_mask = adata.var['highly_variable'].copy()
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_hvg)
        cutoff_mask = adata.var['highly_variable'].copy()
        tf_mask = adata.var_names.isin(all_known_tfs)
        tf_in_default = default_mask & tf_mask
        combined = cutoff_mask.copy()
        if include_hv_tfs:
            combined |= tf_in_default

        if include_hv_tfs:
            adata_subset = SEACell_ad[:, combined].copy()
        else:
            adata_subset = SEACell_ad[:, adata.var['highly_variable']].copy()

        expr_df = pd.DataFrame(adata_subset.X.toarray(), index=adata_subset.obs_names, columns=adata_subset.var_names).T

        filtered_net_df = net_df[
            net_df['Gene1'].isin(adata_subset.var.index) &
            net_df['Gene2'].isin(adata_subset.var.index)
        ].copy()
        filtered_net_df = filtered_net_df[filtered_net_df['Gene1'] != filtered_net_df['Gene2']]
        filtered_net_df.drop_duplicates(keep='first', inplace=True)
        expr_file = param_prefix / 'ExpressionData.csv'
        expr_df.to_csv(expr_file, index=True)
        filtered_network_filename = param_prefix / 'refNetwork.csv'
        filtered_net_df.to_csv(filtered_network_filename, index=False)
        
        ### NAME DEFINED OUTSIDE THE FUNCTION
        datasets.append(param_prefix.name)

        
        end_time = time.time()
        elapsed = end_time - start_time
        print(f"[{datetime.now().strftime('%H:%M:%S')}] END   k={n_SEACells} (took {elapsed:.1f} s)")

        print(f"Saved to {expr_file} and {filtered_network_filename}")
        return {
            'n_cells': adata_subset.n_obs,
            'n_genes': adata_subset.n_vars,
            'k': n_SEACells,
            'model': model
        }




    results_summary = []
    overall_start = time.time()

    results_summary = Parallel(n_jobs=-1)(
        delayed(run_pipeline)(k) for k in k_list
    )

    overall_end = time.time()
    print(f"\n✅ All jobs finished in {overall_end - overall_start:.1f} seconds")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modular BEELINE preprocessing and downsampling script.")
    
    # --- Required Paths ---
    parser.add_argument('--input_data', type=str, required=True, 
                        help="Path to the input .h5ad")
    parser.add_argument('--dataset_name', type=str, required=True, 
                        help="Name for the new dataset (e.g., 'new_data_2025')")
    parser.add_argument('--k_values', type=int, nargs='+', 
                        default = [3, 5, 10, 20, 50, 100, 300],
                        help="List of k values (metacell counts) to run. "
                             "Default: 3, 5, 10, 20, 50, 100, 300")
    # --- Paths with defaults ---
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