#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 







borzoi_effect_file="${1}"
genotype_stem="${2}"
genotype_sample_mapping_file="${3}"
eqtl_sumstats_file="${4}"
anno_method="${5}"
borzoi_based_prior_output_stem="${6}"
calibration_n_bootstraps="${7}"


echo $borzoi_based_prior_output_stem

source ~/.bashrc
conda activate plink_env





if [[ -z "$calibration_n_bootstraps" ]]; then
	python run_borzoi_based_prior_inference_ldscore_grid_based.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file $anno_method $borzoi_based_prior_output_stem
else
	python run_borzoi_based_prior_inference_ldscore_grid_based.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file $anno_method $borzoi_based_prior_output_stem $calibration_n_bootstraps
fi
