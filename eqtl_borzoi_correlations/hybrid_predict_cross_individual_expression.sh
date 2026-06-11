#!/bin/bash
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 






borzoi_effect_file="${1}"
borzoi_annotation_file="${2}"
genotype_stem="${3}"
genotype_sample_mapping_file="${4}"
expr_file="${5}"
distribution="${6}"
borzoi_based_prior_output_stem="${7}"
expr_prediction_output_file="${8}"
standardize_geno="${9}"
training_sample_size="${10}"
max_training_sample_size="${11}"

source ~/.bashrc
conda activate plink_env

echo $expr_prediction_output_file

python hybrid_predict_cross_individual_expression.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno $training_sample_size $max_training_sample_size