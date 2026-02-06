#!/bin/bash
#SBATCH -t 0-5:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                          # Partition to run in
#SBATCH --mem=10GB  



source ~/.bashrc
conda activate borzoi



borzoi_pred_input_stem="${1}"
output_file="${2}"
vg_pairs_to_test_file="${3}"




python organize_borzoi_sed_results.py $borzoi_pred_input_stem $output_file $vg_pairs_to_test_file
