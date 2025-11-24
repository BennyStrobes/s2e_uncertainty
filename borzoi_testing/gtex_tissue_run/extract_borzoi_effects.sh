#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-20:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)



processed_fm_eqtl_output_file="${1}"
borzoi_downloads_dir="${2}"
borzoi_examples_data_dir="${3}"
borzoi_eqtl_output_file="${4}"
sim_iter="${5}"
total_sims="${6}"



source ~/.bashrc
conda activate borzoi


echo "SIMULATION "$sim_iter
date

python extract_borzoi_effects.py $processed_fm_eqtl_output_file $borzoi_downloads_dir $borzoi_examples_data_dir $borzoi_eqtl_output_file $sim_iter $total_sims
date
