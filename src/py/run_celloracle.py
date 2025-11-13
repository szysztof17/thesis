import os
import subprocess
from pathlib import Path
import argparse
import yaml

import celloracle as co
import pandas as pd
import scanpy as sc
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object
import numpy as np

# ==============================================================================
# STEP 1: CICERO ANALYSIS
# ==============================================================================
def run_cicero(config, r_env_path):
    """
    Executes the Cicero R script, including cell filtering and preprocessing.
    """
    print("--- Running Step 1: Cicero Analysis (with preprocessing) ---")
    script_path = Path(config["base_dir"]) / "src/R/cicero_for_celloracle.R"
    output_dir = Path(config["cicero_output_dir"])
    dataset_cicero_dir = output_dir / config["dataset_name"]
    os.makedirs(dataset_cicero_dir, exist_ok=True)

    command = [
        "conda", "run", "-p", r_env_path, "Rscript", str(script_path),
        "--input_dir", str(config["input_atac_dir"]),
        "--dataset_name", str(config["dataset_name"]),
        "--chrom_sizes", str(config["chrom_sizes"]),
        "--output_dir", str(output_dir),
        # Pass the new filtering arguments to the R script
        "--min_count", str(config["min_peaks_per_cell"]),
        "--max_count", str(config["max_peaks_per_cell"])
    ]

    try:
        # Using capture_output to show R script's print statements in real time
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print("Stderr from R script:")
            print(result.stderr)
        print("Cicero analysis completed successfully.")
        return dataset_cicero_dir
    except subprocess.CalledProcessError as e:
        print(f"Error running Cicero R script: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return None

# ==============================================================================
# STEP 2: PROCESS PEAKS
# ==============================================================================
def process_peaks(config, cicero_results_dir):
    """
    Processes Cicero output to get a filtered peak-gene list.
    """
    print("\n--- Running Step 2: Processing Cicero Peaks ---")
    peaks_path = cicero_results_dir / "all_peaks.csv"
    connections_path = cicero_results_dir / "cicero_connections.csv"
    output_path = cicero_results_dir / f"{config['dataset_name']}_peak_gene_list.csv"

    peaks = pd.read_csv(peaks_path, index_col=0).x.values
    cicero_connections = pd.read_csv(connections_path, index_col=0)

    tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=config["ref_genome"])
    integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, cicero_connections=cicero_connections)

    peak_df = integrated[integrated.coaccess >= config["coaccess_threshold"]]
    peak_df = peak_df[["peak_id", "gene_short_name"]].reset_index(drop=True)

    peak_df.to_csv(output_path)
    print(f"Filtered peak-gene list saved to {output_path}")
    return output_path

def preprocess_scrna(config, qc_params):
    """
    Performs QC, filtering, and preparation of scRNA-seq data.
    """
    print("\n--- Running Step 3: Preprocessing scRNA-seq data from raw files ---")
    output_dir = Path(config["beeline_inputs_dir"])
    param_prefix = output_dir / f"{config['dataset_name']}_{config['network_name']}_{config['n_top_hvg']}TFs"
    os.makedirs(param_prefix, exist_ok=True)
    
    adata_path = Path(config["input_adata_path"]) 
    adata = sc.read_10x_h5(adata_path)
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(adata, min_genes=qc_params["min_genes"])
    sc.pp.filter_cells(adata, max_counts=qc_params["max_counts"])
    adata = adata[adata.obs.pct_counts_mt < qc_params["pct_mt"], :].copy()
    adata = adata[:, ~adata.var.mt].copy()

    sc.pp.filter_genes(adata, min_cells=20)
    print(f"Data shape after filtering: {adata.shape}")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=qc_params["n_top_hvg"])
    adata.raw = adata.copy()
    adata.layers['counts']=adata.X.copy() 

    
    sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg)
    
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    adata_hvg_path = param_prefix / "adata.h5ad"
    adata_hvg.write(adata_hvg_path)

    print(f"Processed AnnData object saved to {adata_hvg_path}")

    # --- EXPORT BEELINE-SPECIFIC FILES ---
    genes_in_data = set(adataadataadataadataadata.var_names)
    tfs_df = pd.read_csv(config["tf_list_file"], header=0)
    all_known_tfs = set(tfs_df[tfs_df.columns[0]])
    tfs_in_data = all_known_tfs.intersection(genes_in_data)
    print(f"Found {len(tfs_in_data)} TFs present in the expression data.")

    sc.pp.highly_variable_genes(adataadataadataadata, n_top_genes=None)
    default_mask = adataadataadata.var['highly_variable'].copy()
    n_default = default_mask.sum()
    print(f"Default HVGs: {n_default}")
    # Cutoff HVGs
    sc.pp.highly_variable_genes(adataadata, n_top_genes=n_top_hvg)
    cutoff_mask = adata.var['highly_variable'].copy()
    n_cutoff = cutoff_mask.sum()
    print(f"Top {n_top_hvg} HVGs: {n_cutoff}")
    # Include TFs
    tf_mask = adata.var_names.isin(all_known_tfs)
    tf_in_default = default_mask & tf_mask
    n_tf_default = tf_in_default.sum()
    print(f"TFs in default HVGs: {n_tf_default}")
    combined = cutoff_mask.copy()
    combined |= tf_in_default
    added = combined.sum() - n_cutoff
    print(f"Including additional {added} TFs → total combined HVGs: {combined.sum()}")
    os.makedirs(param_prefix, exist_ok=True)
    # 1. Export ExpressionData.csv
    expr_df = pd.DataFrame(adata_hvg.X.toarray(), index=adata_hvg.obs_names, columns=adata_hvg.var_names).T
    expr_df.to_csv(param_prefix / 'ExpressionData.csv')
    
    # 2. Filter and export Network.csv
    net_df = pd.read_csv(config["network_file"])
    filtered_net_df = net_df[net_df['Gene1'].isin(adata_hvg.var.index) & net_df['Gene2'].isin(adata_hvg.var.index)].copy()
    filtered_net_df.to_csv(param_prefix / 'Network.csv', index=False)

    print(f"BEELINE input files (ExpressionData, Network) saved in: {param_prefix}")


    return adata_hvg_path

def scan_motifs(config, filtered_peaks_path):
    """
    Performs motif scanning on the filtered peaks.
    """
    print("\n--- Running Step 4: Motif Scanning ---")
    output_dir = Path(config["celloracle_output_dir"])
    os.makedirs(output_dir, exist_ok=True)
    genome_installation = ma.is_genome_installed(ref_genome=config["ref_genome"],
                                             genomes_dir=None)
    if not genome_installation:
        import genomepy
        genomepy.install_genome(name=ref_genome, provider="UCSC", genomes_dir=None)
    else:
        print(ref_genome, "is installed.")

    peaks = pd.read_csv(filtered_peaks_path, index_col=0)
    peaks = ma.check_peak_format(peaks, ref_genome=config["ref_genome"])

    tfi = ma.TFinfo(peak_data_frame=peaks, ref_genome=config["ref_genome"])
    tfi.scan(fpr=0.02, motifs=None, verbose=True)

    tfi.reset_filtering()
    tfi.filter_motifs_by_score(threshold=10)
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)
    df = tfi.to_dataframe()
    parquet_path = output_dir / "tf_info.parquet"
    df.to_parquet(parquet_path)

    # If you want, you can save the result as a dictionary as follows.
    td = tfi.to_dictionary(dictionary_type="TF2targetgenes")
    save_as_pickled_object(td, output_dir / "TFinfo_TF2targetgenes.pickled")

    print(f"Motif scan results (tf_info.parquet) saved in: {output_dir}")
    return parquet_path

def run_celloracle(config, adata_path, base_grn):
    """
    Runs the main CellOracle GRN inference using a provided base GRN DataFrame.
    """
    print("\n--- Running Step 5: CellOracle GRN Inference ---")
    output_dir = Path(config["celloracle_output_dir"])
    use_custom_grn = config["use_custom_grn"]
    suffix = 'withATAC' if use_custom_grn else 'noATAC'

    os.makedirs(output_dir, exist_ok=True)
    cluster_column = config["cluster_column"] if config["cluster_column"] else 'leiden'
    adata = sc.read_h5ad(adata_path)
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="leiden", embedding_name="X_umap")
    
    oracle.import_TF_data(TF_info_matrix=base_grn)

    oracle.perform_PCA()

    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    n_comps = min(n_comps, 50)
    n_cell = oracle.adata.shape[0]
    print(f"cell number is :{n_cell}")
    k = int(0.025*n_cell)
    print(f"Auto-selected k is :{k}")
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=16)
    links = oracle.get_links(cluster_name_for_GRN_unit=cluster_column, alpha=10,
                         verbose_level=10)
    dataset_ntopgenes_short = f"{config['dataset_name']}_{config['network_name']}_{config['n_top_hvg']}TFs"

    dataset_ntopgenes = f"{config['dataset_name']}_{config['network_name']}_{config['n_top_hvg']}TFs/{suffix}"
   
    links_path = output_dir / dataset_ntopgenes 
    os.makedirs(links_path, exist_ok=1)
    links_path = str(links_path)
    file_path = links_path+"/links.celloracle.links"
    links.to_hdf5(file_path = file_path)
    links_o = co.load_hdf5(file_path)

    oracle_path = output_dir / dataset_ntopgenes / "oracle.celloracle.oracle"
    oracle_path = str(oracle_path)
    oracle.to_hdf5(str(oracle_path))
    print(f"CellOracle object saved to: {oracle_path}")
    grns_path = output_dir / dataset_ntopgenes  / 'grns'
    os.makedirs(grns_path, exist_ok=1)
    for key in links.links_dict.keys():
        links.links_dict[key].to_csv(f"{grns_path}/raw_GRN_for_{key}.csv")

    directory = Path(grns_path)
    files = [f for f in directory.iterdir() if f.is_file()]

    df = pd.concat([pd.read_csv(f, index_col=0) for f in files], ignore_index=True)
    mask = df.p<0.05
    df = df[mask]
    df = df[~(df.coef_abs==0)] 
    df = df.loc[df.groupby(['source', 'target'])['coef_abs'].idxmax()]
    df = df[['source', 'target', 'coef_abs'] ]
    df.sort_values(by='coef_abs', inplace=True, ascending = False)
    df.reset_index(drop=True, inplace=True)
    df.columns = ['Gene1','Gene2','EdgeWeight']
    df.to_csv(output_dir / dataset_ntopgenes / "rankedEdges.csv", index=False, sep = '\t')


    beeline_outputs_dir = Path(str(config["beeline_inputs_dir"]).replace('inputs', 'outputs'))
    os.makedirs(beeline_outputs_dir, exist_ok=1)
    method_name= 'CellOracle_withATAC' if use_custom_grn else 'CellOracle_noATAC'
    beeline_oupputs_method_path = beeline_outputs_dir / dataset_ntopgenes_short / method_name
    os.makedirs(beeline_oupputs_method_path, exist_ok = True)
    df.to_csv(beeline_oupputs_method_path / "rankedEdges.csv", index=False, sep = '\t')

# ==============================================================================
# MAIN ORCHESTRATION BLOCK
# ==============================================================================
def main():
    # --- 1. Set up argument parser to read config file path ---
    parser = argparse.ArgumentParser(description="Run the CellOracle pipeline from a configuration file.")
    parser.add_argument("config", help="Path to the YAML configuration file.")
    args = parser.parse_args()

    # --- 2. Load and process the configuration file ---
    with open(args.config, 'r') as file:
        user_config = yaml.safe_load(file)

    # Convert string paths to Path objects and build the final config dictionary
    base_dir = Path(user_config["BASE_DIR"])

    config = {
        "base_dir": base_dir,
        "dataset_name": user_config["DATASET_NAME"],
        'network_name': user_config["RAW_INPUT_DIRS"]["NETWORK_NAME"],
        'network_file': user_config["RAW_INPUT_DIRS"]["NETWORK_FILE"],
        "tf_list_file": base_dir / user_config["RAW_INPUT_DIRS"]["TF_LIST"],
        "ref_genome": user_config["REF_GENOME"],
        # Paths derived from base_dir and dataset_name
        "input_atac_dir": base_dir / user_config["RAW_INPUT_DIRS"]["ATAC"],
        "input_rna_dir": base_dir / user_config["RAW_INPUT_DIRS"]["RNA"],
        "input_adata_path": base_dir / user_config["RAW_INPUT_DIRS"]["H5"],
        "chrom_sizes": base_dir / user_config["RAW_INPUT_DIRS"]["CHROM_SIZES"],
        "cluster_column": user_config["CLUSTER_COLUMN"],
        "cicero_output_dir": base_dir / "results/cicero_output",
        "celloracle_output_dir": base_dir / "results/celloracle_output" / user_config["DATASET_NAME"],
        "beeline_inputs_dir": base_dir / "BEELINE/inputs" / user_config["DATASET_NAME"],
        # QC and Threshold parameters
        "min_peaks_per_cell": user_config["MIN_PEAKS_PER_CELL"],
        "max_peaks_per_cell": user_config["MAX_PEAKS_PER_CELL"],
        "coaccess_threshold": user_config["COACCESS_THRESHOLD"],
        "min_genes_per_cell": user_config["MIN_GENES_PER_CELL"],
        "max_counts_per_cell": user_config["MAX_COUNTS_PER_CELL"],
        "pct_mt_per_cell": user_config["PCT_MT_PER_CELL"],
        "n_top_hvg": user_config["N_TOP_HVG"],
        "use_custom_grn": user_config["USE_ATAC_SEQ"] 
    }

    # --- 3. Run Pipeline based on configuration ---
    base_grn_to_use = None
    if user_config["USE_ATAC_SEQ"]:
        print("PIPELINE MODE: Using custom GRN from ATAC-seq data.")
        
        # New logic to check for intermediate files
        custom_grn_path = user_config.get("ATAC_INTERMEDIATE_PATHS", {}).get("CUSTOM_GRN")
        peak_list_path = user_config.get("ATAC_INTERMEDIATE_PATHS", {}).get("PEAK_GENE_LIST")

        if custom_grn_path:
            print(f"--- SKIPPING ATAC processing: Loading pre-computed GRN from {custom_grn_path} ---")
            base_grn_to_use = pd.read_parquet(base_dir / custom_grn_path)
        elif peak_list_path:
            print(f"--- SKIPPING Cicero & Peak Processing: Using peak list from {peak_list_path} ---")
            base_grn_to_use = scan_motifs(config, base_dir / peak_list_path)
        else:
            print("--- Running full ATAC-seq workflow from raw data ---")
            cicero_dir = run_cicero(config, user_config["R_CONDA_ENV_PATH"])
            if cicero_dir:
                peaks_file = process_peaks(config, cicero_dir)
                base_grn_to_use = scan_motifs(config, peaks_file)
    else:
        print("PIPELINE MODE: Using default promoter-based GRN.")
        if config["ref_genome"] == "hg38":
            base_grn_to_use = co.data.load_human_promoter_base_GRN(version = 'hg38_gimmemotifsv5_fpr2')
        elif config["ref_genome"] == "hg19":
            base_grn_to_use = co.data.load_human_promoter_base_GRN(version = 'hg19_gimmemotifsv5_fpr2')
        else:
            raise ValueError("Default GRN only available for 'hg38' or 'hg19'.")
    if base_grn_to_use is not None:
        adata_path_str = user_config.get("PREPROCESSED_ADATA_PATH")
        if adata_path_str:
            adata_file = base_dir / adata_path_str
            print(f"\nUsing pre-processed AnnData from {adata_file}")
        else:
            adata_file = preprocess_scrna(config)
        
        run_celloracle(config, adata_file, base_grn_to_use)
        print("\n🎉 Pipeline finished successfully! 🎉")
    else:
        print("\nPipeline stopped because a base GRN could not be generated or loaded.")

if __name__ == "__main__":
    main()






