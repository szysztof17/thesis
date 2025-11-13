#!/usr/bin/env Rscript
library(cicero)
library(monocle3)
library(Matrix)
library(argparse)

parser <- ArgumentParser(description='Run Cicero analysis')
parser$add_argument('--input_dir', help='Directory with matrix.mtx.gz, barcodes.tsv.gz, and peaks.bed.gz files')
parser$add_argument('--dataset_name', help='Name of the dataset')
parser$add_argument('--chrom_sizes', help='Path to chromosome sizes file')
parser$add_argument('--output_dir', help='Output directory for Cicero results')

parser$add_argument('--min_count', type="integer", default=1000, help='Minimum total counts per cell')
parser$add_argument('--max_count', type="integer", default=30000, help='Maximum total counts per cell')

args <- parser$parse_args()

# --- Configuration ---
dataset_name <- args$dataset_name
data_dir <- args$input_dir
chrom_sizes_file <- args$chrom_sizes
output_folder <- file.path(args$output_dir, dataset_name)
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# --- Data Loading ---
in_mtx <- file.path(data_dir, "matrix.mtx.gz")
in_barcodes <- file.path(data_dir, "barcodes.tsv.gz")
in_bed <- file.path(data_dir, "peaks.bed.gz")

indata <- Matrix::readMM(in_mtx)
indata@x[indata@x > 0] <- 1

cellinfo <- read.table(in_barcodes, header = FALSE)
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

peakinfo <- read.table(in_bed)
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# --- Create CDS object ---
input_cds <- suppressWarnings(new_cell_data_set(indata,
                                                cell_metadata = cellinfo,
                                                gene_metadata = peakinfo))
input_cds <- monocle3::detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]


# 1. Filter cells based on library size (total counts)
print(paste("Original number of cells:", ncol(input_cds)))
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= args$min_count]
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= args$max_count]
print(paste("Number of cells after filtering:", ncol(input_cds)))

# --- Run Cicero ---
set.seed(420)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

conns <- run_cicero(cicero_cds, chrom_sizes_file)

# --- Save results ---
write.csv(x = conns, file = file.path(output_folder, "cicero_connections.csv"))
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = file.path(output_folder, "all_peaks.csv"))

print("Cicero analysis complete.")