#!/usr/bin/env python3
import os
# Disable Numba's Just-In-Time compiler
os.environ['NUMBA_DISABLE_JIT'] = '1'

# Force single-threaded operation to avoid thread-counting crashes
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['NUMBA_NUM_THREADS'] = '1'
#os.environ['NUMBA_THREADING_LAYER'] = 'workqueue'


import os
import yaml
import subprocess
import pickle
import pandas as pd
import scanpy as sc
from pathlib import Path
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl
from pycisTopic.lda_models import run_cgs_models_mallet
from pycisTopic.lda_models import evaluate_models
from pycisTopic.topic_binarization import *
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
from pycisTopic.diff_features import *




MACS_PATH = "/home/kl467102/miniconda3/envs/scenic_plus/bin/macs2"
PYCISTOPIC = "/home/kl467102/miniconda3/envs/scenic_plus/bin/pycistopic"
MALLET = "Mallet-202108/bin/mallet"

def log(msg):
    print(f"[SCENIC+] {msg}")

def run_pycistopic(cfg):
    """Stage 1: pycisTopic preprocessing and consensus peaks"""
    variable_name = cfg["inputs"]["variable_name"]

    ds = cfg["dataset"]
    outdir = Path(cfg["outdir"])
    tmpdir = outdir / "tmp"
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)

    log(f"Running pycisTopic preprocessing for {ds}...")
    adata = sc.read_h5ad(cfg["inputs"]["rna"])
    cell_data = adata.obs.copy()
    cell_data["label"] = (
        cell_data["label"]
        .str.replace(" ", "_")
        .str.replace("(", "")
        .str.replace(")", "")
    )
    cell_data["sample_id"] = ds

    chromsizes = pd.read_csv(
        cfg["inputs"]["chromsizes"],
        sep="\t", header=None, names=["Chromosome", "End"]
    )
    chromsizes.insert(1, "Start", 0)

    fragments_dict = {ds: cfg["inputs"]["fragments"]}

    bed_dir = outdir / "beds"
    bw_dir = outdir / "bws"
    os.makedirs(bed_dir, exist_ok=True)
    os.makedirs(bw_dir, exist_ok=True)

    bw_paths, bed_paths = export_pseudobulk(
        input_data=cell_data,
        variable=variable_name,
        sample_id_col="sample_id",
        chromsizes=chromsizes,
        bed_path=str(bed_dir),
        bigwig_path=str(bw_dir),
        path_to_fragments=fragments_dict,
        n_cpu=cfg["n_cpu"],
        normalize_bigwig=True,
        temp_dir=str(tmpdir),
        split_pattern = "-"
    )

    log("Peak calling...")
    macs_output_dir = outdir / "MACS"
    narrow_peaks_dict = peak_calling(
        macs_path=MACS_PATH,
        bed_paths=bed_paths,
        outdir=str(macs_output_dir),
        genome_size="hs",
        n_cpu=cfg["n_cpu"],
    )
    path_to_blacklist = cfg["inputs"]["blacklist"]
    log("Generating consensus peaks...")
    consensus_peaks = get_consensus_peaks(
        narrow_peaks_dict=narrow_peaks_dict,
        peak_half_width=250,
        chromsizes=chromsizes,
        path_to_blacklist=cfg["inputs"]["blacklist"],
    )

    consensus_bed = outdir / "consensus_regions.bed"
    consensus_peaks.to_bed(path=consensus_bed, keep=True)

    path_to_regions = {ds : consensus_bed}

    pycistopic_qc_output_dir = outdir / 'qc'
    os.makedirs(pycistopic_qc_output_dir, exist_ok=True)

    cmd = (
        f"{PYCISTOPIC} qc "
        f" --fragments {fragments_dict[ds]} "
        f" --regions {path_to_regions[ds]} "
        f' --tss {cfg["inputs"]["tss"]} '
        f" --output {pycistopic_qc_output_dir}/{ds}"
    )

    result = subprocess.run(
        cmd,
        shell=True,           # needed for the shell to interpret the command correctly
        capture_output=True,  # capture stdout and stderr
        text=True             # return strings instead of bytes
    )

    # Print standard output and standard error
    print(result.stdout)
    print(result.stderr)

    sample_id_to_barcodes_passing_filters = {}
    sample_id_to_thresholds = {}
    for sample_id in fragments_dict:
        (
            sample_id_to_barcodes_passing_filters[sample_id],
            sample_id_to_thresholds[sample_id]
        ) = get_barcodes_passing_qc_for_sample(
                sample_id = sample_id,
                pycistopic_qc_output_dir = pycistopic_qc_output_dir,
                unique_fragments_threshold = None, # use automatic thresholding
                tss_enrichment_threshold = None, # use automatic thresholding
                frip_threshold = 0,
                use_automatic_thresholds = True,
        )



    path_to_regions = path_to_regions[ds]

    cistopic_obj_list = []
    for sample_id in fragments_dict:
        sample_metrics = pl.read_parquet(
            os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
        ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments = fragments_dict[sample_id],
            path_to_regions = str(path_to_regions),
            path_to_blacklist = path_to_blacklist,
            metrics = sample_metrics,
            valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
            n_cpu = cfg["n_cpu"],
            project = sample_id,
            split_pattern = '-'
        )
        cistopic_obj_list.append(cistopic_obj)
    cistopic_obj = cistopic_obj_list[0]


    cistopic_obj.add_cell_data(cell_data, split_pattern='-')
    pickle.dump(cistopic_obj, open(outdir / "cistopic_obj.pkl", "wb"))
    log(f"Saved cistopic object at {outdir}/cistopic_obj.pkl")

    os.environ['MALLET_MEMORY'] = '200G'

    save_path = outdir / "mallet_models"
    os.makedirs(save_path, exist_ok=True)
    models=run_cgs_models_mallet(
        cistopic_obj,
        #n_topics=[2, 5, 10, 15, 30,  40, 50],
        n_topics=[ 40],

        n_cpu=cfg["n_cpu"],
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        save_path=save_path,
        mallet_path=MALLET,
    )

    model = evaluate_models(
        models,
        select_model = 40,
        return_model = True
    )
    
    cistopic_obj.add_LDA_model(model)


    from pycisTopic.clust_vis import run_umap
    run_umap(cistopic_obj, target  = 'cell', scale=True, verbose=False)

    region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
    region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)
    binarized_cell_topic = binarize_topics(
        cistopic_obj,
        target='cell',
        method='li')

    topic_annot = topic_annotation(
        cistopic_obj,
        annot_var='label',
        binarized_cell_topic=binarized_cell_topic,
        general_topic_thr=0.2
    )

    imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
    variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
    markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='label', var_features=variable_regions, split_pattern = '-', n_cpu=24)

    candidate_enhancers_path = outdir / 'candidate_enhancers'
    os.makedirs(candidate_enhancers_path, exist_ok=True)
    folder = os.path.join(outdir, 'region_sets', 'Topics_otsu')
    os.makedirs(folder, exist_ok=True)
    for topic in binarized_cell_topic:
        region_names_to_coordinates(
            region_bin_topics_otsu[topic].index
        ).sort_values(
            ['Chromosome', 'Start', 'End']
        ).to_csv(
            os.path.join(folder, f'{topic}.bed'),
            sep='\t',
            header=False, index=False
        )
    folder = os.path.join(outdir, 'region_sets', 'Topics_top_3k')
    os.makedirs(folder, exist_ok=True)
    for topic in region_bin_topics_top3k:
        region_names_to_coordinates(
            region_bin_topics_top3k[topic].index
        ).sort_values(
            ['Chromosome', 'Start', 'End']
        ).to_csv(
            os.path.join(folder, f'{topic}.bed'),
            sep='\t',
            header=False, index=False
        )

    # Save DARs
    folder = os.path.join(outdir, 'region_sets', 'DARs_cell_type')
    os.makedirs(folder, exist_ok=True)
    for cell_type in markers_dict:
        region_names_to_coordinates(
            markers_dict[cell_type].index
        ).sort_values(
            ['Chromosome', 'Start', 'End']
        ).to_csv(
            os.path.join(folder, f'{cell_type}.bed'),
            sep='\t',
            header=False,
            index=False
        )

def run_snakemake(cfg):
    """Stage 2: Run SCENIC+ Snakemake pipeline"""
    snk = cfg["snakemake"]
    log("Launching Snakemake pipeline...")

    cmd = [
        "snakemake",
        "-s", snk["snakefile"],
        "--configfile", snk["configfile"],
        "--cores", str(snk["cores"]),
    ]
    if snk.get("use_conda", True):
        cmd.append("--use-conda")
    if snk.get("rerun_incomplete", True):
        cmd.append("--rerun-incomplete")

    log(" ".join(cmd))
    subprocess.run(cmd, check=True)
    log("Snakemake completed successfully.")

def postprocess(cfg):
    """Stage 3: Extract ranked edges from SCENIC+ output"""
    import mudata
    log("Extracting GRN edges from SCENIC+ results...")

    mdata_path = Path(cfg["outdir"]) / "Snakemake" / "outs" / "scplusmdata.h5mu"
    if not mdata_path.exists():
        raise FileNotFoundError(f"Cannot find {mdata_path}")

    mdata = mudata.read(mdata_path)
    outdir = Path(cfg["outdir"]) / "results"
    outdir.mkdir(exist_ok=True, parents=True)

    eReg = pd.concat(
        [mdata.uns["direct_e_regulon_metadata"], mdata.uns["extended_e_regulon_metadata"]],
        axis=0,
    )
    eReg = eReg[['TF', 'Gene', 'rho_TF2G', 'Region']].drop_duplicates()
    eReg = eReg[eReg.rho_TF2G != 0]
    eReg.columns = ["Gene1", "Gene2", "EdgeWeight", "peak"]
    net = eReg.groupby(["Gene1", "Gene2"], as_index=False)["EdgeWeight"].max()
    net["EdgeWeight"] = net["EdgeWeight"].abs()
    net.sort_values("EdgeWeight", ascending=False, inplace=True)
    net.to_csv(outdir / "rankedEdges.csv", sep="\t", index=False)
    log(f"Exported rankedEdges.csv to {outdir}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Unified SCENIC+ wrapper")
    parser.add_argument("--config", required=True, help="YAML configuration file")
    parser.add_argument("--run-all", action="store_true", help="Run all stages")
    parser.add_argument("--only", choices=["pycistopic", "snakemake", "postprocess"], help="Run only a specific stage")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    if args.run_all or (args.only == "pycistopic"):
        if cfg["options"].get("run_pycistopic", True):
            run_pycistopic(cfg)
    if args.run_all or (args.only == "snakemake"):
        if cfg["options"].get("run_snakemake", True):
            run_snakemake(cfg)
    if args.run_all or (args.only == "postprocess"):
        if cfg["options"].get("postprocess", True):
            postprocess(cfg)

if __name__ == "__main__":
    main()
