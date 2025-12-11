#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-9:30                         # Runtime in D-HH:MM format
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)

source ~/.bashrc
conda activate borzoi



raw_fine_mapping_eqtl_results_file="$1"
pip_threshold="$2"
processed_fm_eqtl_output_file="$3"
fm_vcf_output_file="${4}"

python parse_eqtl_data.py $raw_fine_mapping_eqtl_results_file $pip_threshold $processed_fm_eqtl_output_file

python make_vcf_out_of_eqtls.py $processed_fm_eqtl_output_file $fm_vcf_output_file

echo $fm_vcf_output_file