#!/bin/bash
#SBATCH -t 0-34:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=30GB 







borzoi_effect_file="${1}"
borzoi_annotation_file="${2}"
genotype_stem="${3}"
genotype_sample_mapping_file="${4}"
expr_file="${5}"
distribution="${6}"
borzoi_based_prior_output_stem="${7}"
standardize_geno="${8}"
maxiter="${9}"


echo $borzoi_based_prior_output_stem

source ~/.bashrc
conda activate plink_env


if [[ "$distribution" == "gaussian" ]]; then
	python run_borzoi_based_prior_inference.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
elif [[ "$distribution" == "ashr" ]]; then
	python run_borzoi_based_prior_inference_ashr.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
elif [[ "$distribution" == "effect_size_grid" ]]; then
	python run_borzoi_based_prior_inference_effect_size_grid.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem
elif [[ "$distribution" == "uniform_prior_grid" ]]; then
	python run_borzoi_based_prior_inference_uniform_grid.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
elif [[ "$distribution" == "glm_borzoi_magnitude" ]]; then
	python run_borzoi_based_prior_inference_glm_borzoi_magnitude.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
elif [[ "$distribution" == "borzoi_scaled_variance" ]]; then
	python run_borzoi_based_prior_inference_borzoi_scaled_variance.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
elif [[ "$distribution" == "glm_borzoi_magnitude_marginalized" ]]; then
	if [[ -z "$maxiter" ]]; then
		python run_borzoi_based_prior_inference_glm_borzoi_magnitude_marginalized.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
	else
		python run_borzoi_based_prior_inference_glm_borzoi_magnitude_marginalized.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno $maxiter
	fi


else
	echo "Unknown distribution: $distribution"
	exit 1
fi
