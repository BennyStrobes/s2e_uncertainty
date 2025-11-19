#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-1:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)

source ~/.bashrc
conda activate SAGEnet



raw_fine_mapping_eqtl_results_file="$1"
tissue_name="$2"
pip_threshold="$3"
processed_fm_eqtl_output_file="$4"

python parse_eqtl_data.py $raw_fine_mapping_eqtl_results_file $tissue_name $pip_threshold $processed_fm_eqtl_output_file