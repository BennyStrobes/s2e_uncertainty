#!/bin/bash
#SBATCH -t 0-14:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=4GB  



borzoi_gtex_predictions="${1}"
full_gtex_target_file="${2}"
organized_borzoi_gtex_predictions="${3}"
processed_fm_eqtl_output_file="${4}"
gene_tss_file="${5}"


source ~/.bashrc
conda activate borzoi


if false; then
python organize_gtex_borzoi_sed_results_fm.py $borzoi_gtex_predictions $full_gtex_target_file $organized_borzoi_gtex_predictions
fi



if false; then
python organize_gtex_borzoi_sed_results_fm_gene_span_debug.py $borzoi_gtex_predictions $full_gtex_target_file $organized_borzoi_gtex_predictions
fi

if false; then
python organize_gtex_borzoi_sed_results_fm_random_rc_debug.py $borzoi_gtex_predictions $full_gtex_target_file $organized_borzoi_gtex_predictions
fi



python extract_fine_mapped_effects.py $full_gtex_target_file $processed_fm_eqtl_output_file $organized_borzoi_gtex_predictions $gene_tss_file
