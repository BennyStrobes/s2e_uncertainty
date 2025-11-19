#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-1:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)

source ~/.bashrc
conda activate SAGEnet




gtex_wb_tpm_file="${1}"
gtex_wb_processed_expression_file_stem="${2}"
tpm_thresh="${3}"
max_prop_sample_missing="${4}"

python process_gtex_expression_data.py $gtex_wb_tpm_file $gtex_wb_processed_expression_file_stem $tpm_thresh $max_prop_sample_missing