#!/bin/bash
#SBATCH -t 0-2:00
#SBATCH -p bch-compute
#SBATCH --mem=20GB

gtex_fm_file="${1}"
target_tissue="${2}"
fine_mapping_method="${3}"
pip_threshold="${4}"
borzoi_effect_file="${5}"
borzoi_based_prior_output_stem="${6}"
output_file="${7}"

source ~/.bashrc
conda activate plink_env

python compare_finemapped_eqtls_to_borzoi_prior.py $gtex_fm_file $target_tissue $fine_mapping_method $pip_threshold $borzoi_effect_file $borzoi_based_prior_output_stem $output_file
