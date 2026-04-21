#!/bin/bash
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 



input_expression_file="${1}"
input_covariate_file="${2}"
residual_expression_file="${3}"
gtex_v10_pc_genes_gtf="${4}"



source ~/.bashrc
conda activate plink_env


python residualize_expression.py $input_expression_file $input_covariate_file $residual_expression_file $gtex_v10_pc_genes_gtf