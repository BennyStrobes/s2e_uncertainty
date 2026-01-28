#!/bin/bash
#SBATCH -t 0-60:50                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=10GB



sig_thresh="${1}"
tissue_overlap_output_file="${2}"
gtex_tissue_names_file="${3}"
organized_borzoi_gtex_predictions="${4}"


source ~/.bashrc
conda activate borzoi


python compute_tissue_overlap_statistics.py $sig_thresh $tissue_overlap_output_file $gtex_tissue_names_file $organized_borzoi_gtex_predictions