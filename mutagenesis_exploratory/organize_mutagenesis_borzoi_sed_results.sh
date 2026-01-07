#!/bin/bash
#SBATCH -t 0-75:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=4GB  



borzoi_gtex_predictions="${1}"
full_gtex_target_file="${2}"
organized_borzoi_gtex_predictions="${3}"
processed_fm_eqtl_output_file="${4}"

source ~/.bashrc
conda activate borzoi

if false; then
python organize_gtex_borzoi_sed_results_fm.py $borzoi_gtex_predictions $full_gtex_target_file $organized_borzoi_gtex_predictions
fi

python extract_mutagenesis_effects.py $full_gtex_target_file $processed_fm_eqtl_output_file $organized_borzoi_gtex_predictions
