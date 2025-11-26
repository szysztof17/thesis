import argparse
import os
import pandas as pd
import scipy
import scanpy as sc
import time

from LingerGRN.preprocess import *
from LingerGRN.pseudo_bulk import *
import LingerGRN.LINGER_tr as LINGER_tr
import LingerGRN.LL_net as LL_net

# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Run LINGER analysis on RNA and ATAC data.")
    
    # Input and output directories
    parser.add_argument('--input_dir', type=str, required=True, help="Path to 10x filtered feature barcodes matrix directory")
    parser.add_argument('--label_file', type=str, help="Path to cell-type annotations file (default: 'label.txt' in the input directory)")
    parser.add_argument('--Datadir', type=str, required=True, help="Directory for the downloaded general gene regulatory network")
    parser.add_argument('--genome', type=str, required=True, choices=['hg38', 'mm10'], help="Genome version (hg38 or mm10)")

    # Optional arguments
    parser.add_argument('--method', type=str, default='LINGER', help="Method for GRN analysis (default: LINGER)")
    parser.add_argument('--outdir', type=str, default='default', help="Directory for the output results (default: current_directory/results/)")
    parser.add_argument('--linger_storage', type=str, default='default', help="Directory to store preprocessed data")

    args = parser.parse_args()

    # Set default label file if not provided
    if not args.label_file:
        args.label_file = os.path.join(args.input_dir, 'label.txt')
    if args.outdir == 'default':
        args.outdir = os.path.join(os.getcwd(),'output/')
    if args.linger_storage == 'default':
        args.linger_storage = os.path.join(os.getcwd(),'data/')
    return args

# Function to track time
def log_time(start_time, step_name):
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"{step_name} completed in {elapsed_time:.2f} seconds.")
    return end_time

# Main function to perform the LINGER analysis
def main():
    # Parse the command-line arguments
    args = parse_args()

    # Initialize start time for the entire process
    start_time = time.time()

    # Step 1 - Preprocessing
    print("\nStep 1: Preprocessing data...")
    step_start_time = time.time()

    input_dir = args.input_dir
    label_file = args.label_file
    linger_storage = args.linger_storage
    Datadir = args.Datadir
    GRNdir = Datadir + 'data_bulk/'

    genome = args.genome
    outdir = args.outdir
    method = args.method
    os.makedirs(outdir, exist_ok = True) 

    # Read in the data
    matrix_path = os.path.join(input_dir, 'matrix.mtx')

    matrix = scipy.io.mmread(matrix_path)
    features = pd.read_csv(os.path.join(input_dir, 'features.tsv'), sep='\t', header=None)
    barcodes = pd.read_csv(os.path.join(input_dir, 'barcodes.tsv'), sep='\t', header=None)
    label = pd.read_csv(label_file, sep='\t', header=0)
    print("Data read successfully.")
    print(f"{pd.concat([features.head(2), features.tail(2)])}")
    step_end_time = log_time(step_start_time, "Data reading")

    step_start_time = time.time()

    # Extract matrices for RNASeq and ATACSeq
    adata_RNA, adata_ATAC = get_adata(matrix, features, barcodes, label)

    # Remove low count cells
    sc.pp.filter_cells(adata_RNA, min_genes=200)
    sc.pp.filter_genes(adata_RNA, min_cells=3)
    sc.pp.filter_cells(adata_ATAC, min_genes=200)
    sc.pp.filter_genes(adata_ATAC, min_cells=3)

    # Find common barcodes between RNA and ATAC datasets
    selected_barcode = list(set(adata_RNA.obs['barcode']) & set(adata_ATAC.obs['barcode']))
    adata_RNA = adata_RNA[adata_RNA.obs['barcode'].isin(selected_barcode)]
    adata_ATAC = adata_ATAC[adata_ATAC.obs['barcode'].isin(selected_barcode)]

    # Generate metacells
    samplelist = list(set(adata_ATAC.obs['sample'].values))
    TG_pseudobulk = pd.DataFrame([])
    RE_pseudobulk = pd.DataFrame([])

    singlepseudobulk = (adata_RNA.obs['sample'].unique().shape[0]*adata_RNA.obs['sample'].unique().shape[0] > 100)
    for tempsample in samplelist:
        adata_RNAtemp = adata_RNA[adata_RNA.obs['sample'] == tempsample]
        adata_ATACtemp = adata_ATAC[adata_ATAC.obs['sample'] == tempsample]
        TG_pseudobulk_temp, RE_pseudobulk_temp = pseudo_bulk(adata_RNAtemp, adata_ATACtemp, singlepseudobulk)
        TG_pseudobulk = pd.concat([TG_pseudobulk, TG_pseudobulk_temp], axis=1)
        RE_pseudobulk = pd.concat([RE_pseudobulk, RE_pseudobulk_temp], axis=1)
        RE_pseudobulk[RE_pseudobulk > 100] = 100

    # Write preprocessed data to files
    os.makedirs(linger_storage, exist_ok=True)
    adata_ATAC.write(os.path.join(linger_storage, 'adata_ATAC.h5ad'))
    adata_RNA.write(os.path.join(linger_storage, 'adata_RNA.h5ad'))
    TG_pseudobulk = TG_pseudobulk.fillna(0)
    RE_pseudobulk = RE_pseudobulk.fillna(0)
    os.makedirs('data', exist_ok=True)

    pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv('data/Peaks.txt', header=None, index=None)
    TG_pseudobulk.to_csv(os.path.join(linger_storage, 'TG_pseudobulk.tsv'))
    RE_pseudobulk.to_csv(os.path.join(linger_storage, 'RE_pseudobulk.tsv'))

    step_end_time = log_time(step_start_time, "Preprocessing")

    # Step 2 - Training
    print("\nStep 2: Training model...")
    step_start_time = time.time()

    preprocess(TG_pseudobulk, RE_pseudobulk, GRNdir, genome, method, outdir)
    activef = 'ReLU'
    LINGER_tr.training(GRNdir, method, outdir, activef, 'Human')

    step_end_time = log_time(step_start_time, "Training model")

    # Step 3 - Population GNR Inference
    print("\nStep 3: Population GNR inference...")
    step_start_time = time.time()

    LL_net.TF_RE_binding(GRNdir, adata_RNA, adata_ATAC, genome, method, outdir)
    LL_net.cis_reg(GRNdir, adata_RNA, adata_ATAC, genome, method, outdir)
    LL_net.trans_reg(GRNdir, method, outdir, genome)

    step_end_time = log_time(step_start_time, "Population GNR inference")

    # Step 4 - Cell-type Specific GNR
    print("\nStep 4: Cell-type specific GNR inference...")
    step_start_time = time.time()

    celltype = 'all'  # Use 'all' or specify a particular label
    LL_net.cell_type_specific_TF_RE_binding(GRNdir, adata_RNA, adata_ATAC, genome, celltype, outdir, method)
    LL_net.cell_type_specific_cis_reg(GRNdir, adata_RNA, adata_ATAC, genome, celltype, outdir, method)
    LL_net.cell_type_specific_trans_reg(GRNdir, adata_RNA, celltype, outdir)

    step_end_time = log_time(step_start_time, "Cell-type specific GNR inference")

    # Final message with total elapsed time
    total_time = log_time(start_time, "LINGER analysis")
    print(f"Total time taken for the analysis: {total_time:.2f} seconds.")

if __name__ == "__main__":
    main()
