#!/bin/bash
#SBATCH -t 0-30:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=30GB 



borzoi_gtex_independent_target_names_file="${1}"
gtex_v10_pc_genes_gtf="${2}"
eqtl_sumstats_dir="${3}"
borzoi_results_dir="${4}"
genotype_stem="${5}"
bayes_input_data_stem="${6}"

source ~/.bashrc
conda activate plink_env


python preprocess_data_for_bayesian_ld_corr.py $borzoi_gtex_independent_target_names_file $gtex_v10_pc_genes_gtf $eqtl_sumstats_dir $borzoi_results_dir $genotype_stem $bayes_input_data_stem