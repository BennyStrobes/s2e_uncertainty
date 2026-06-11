#!/bin/bash
#SBATCH -t 0-4:00
#SBATCH -p bch-compute
#SBATCH --mem=30GB

borzoi_effect_file="${1}"
genotype_stem="${2}"
genotype_sample_mapping_file="${3}"
expr_file="${4}"
borzoi_based_prior_output_stem="${5}"
output_file="${6}"
ld_pruned_output_file="${7}"
ld_prune_r_sq_threshold="${8}"
ld_prune_random_seed="${9}"
distribution="${10:-auto}"

source ~/.bashrc
conda activate plink_env

python predict_marginal_eqtl_effects_from_borzoi_prior.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $borzoi_based_prior_output_stem $output_file $ld_pruned_output_file $ld_prune_r_sq_threshold $ld_prune_random_seed $distribution
