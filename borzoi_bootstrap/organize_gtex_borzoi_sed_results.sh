#!/bin/bash
#SBATCH -t 0-18:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=3GB  



borzoi_gtex_predictions="${1}"
full_gtex_target_file="${2}"
organized_borzoi_gtex_predictions="${3}"

source ~/.bashrc
conda activate borzoi

python organize_gtex_borzoi_sed_results.py $borzoi_gtex_predictions $full_gtex_target_file $organized_borzoi_gtex_predictions