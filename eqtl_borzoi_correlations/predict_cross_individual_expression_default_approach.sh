#!/bin/bash
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 




borzoi_effect_file="${1}"
genotype_stem="${2}"
genotype_sample_mapping_file="${3}"
expr_file="${4}"
expr_prediction_output_file="${5}"
gene_individual_prediction_output_file="${6}"

source ~/.bashrc
conda activate plink_env


python predict_cross_individual_expression_default_approach.py $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $expr_prediction_output_file $gene_individual_prediction_output_file



