#!/bin/bash
#SBATCH -t 0-40:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 






bayes_input_data_summary_file="${1}"
tissue_name="${2}"
borzoi_gtex_independent_target_names_file="${3}"
bayes_ld_corr_output_stem="${4}"

source ~/.bashrc
conda activate plink_env


echo $bayes_ld_corr_output_stem

python run_bayes_ld_corr_magnitude_stratified.py $bayes_input_data_summary_file $tissue_name $borzoi_gtex_independent_target_names_file $bayes_ld_corr_output_stem
