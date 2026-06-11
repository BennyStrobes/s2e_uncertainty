#!/bin/bash
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 







borzoi_effect_file="${1}"
genotype_stem="${2}"
genotype_sample_mapping_file="${3}"
eqtl_sumstats_file="${4}"
anno_method="${5}"
borzoi_based_prior_output_stem="${6}"
model_effect_scale="${7}"
calibration_n_bootstraps="${8}"


echo $borzoi_based_prior_output_stem

source ~/.bashrc
conda activate plink_env





if [[ -z "$model_effect_scale" ]]; then
	python run_borzoi_based_prior_inference_ldscore_based.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file $anno_method $borzoi_based_prior_output_stem
elif [[ -z "$calibration_n_bootstraps" ]]; then
	python run_borzoi_based_prior_inference_ldscore_based.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file $anno_method $borzoi_based_prior_output_stem $model_effect_scale
else
	python run_borzoi_based_prior_inference_ldscore_based.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file $anno_method $borzoi_based_prior_output_stem $model_effect_scale $calibration_n_bootstraps
fi
