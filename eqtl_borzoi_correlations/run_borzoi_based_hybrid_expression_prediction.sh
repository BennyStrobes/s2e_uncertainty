#!/bin/bash
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
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
fold_number="${9}"
n_folds="${10}"
down_sampling_fraction="${11}"

date
echo $borzoi_based_prior_output_stem

source ~/.bashrc
conda activate plink_env



if [[ "$distribution" == "gaussian" ]]; then
	python run_borzoi_based_hybrid_expression_prediction.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno $fold_number $n_folds $down_sampling_fraction
elif [[ "$distribution" == "elastic_net" ]]; then
	python run_elastic_net_expression_prediction.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno $fold_number $n_folds $down_sampling_fraction
elif [[ "$distribution" == "no_borzoi_genome_wide_ridge" ]]; then
	python run_no_borzoi_genome_wide_ridge_expression_prediction.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno $fold_number $n_folds $down_sampling_fraction


else
	echo "Unknown distribution: $distribution"
	exit 1
fi


date