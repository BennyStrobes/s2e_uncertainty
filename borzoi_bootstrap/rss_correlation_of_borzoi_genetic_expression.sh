#!/bin/bash
#SBATCH -t 0-18:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=3GB  





tissue_name="${1}"
borzoi_effect_size_file="${2}"
eqtl_sumstats_dir="${3}"
protein_coding_genes_file="${4}"
genotype_1000G_plink_stem="${5}"
hm3_snp_list_file="${6}"


source ~/.bashrc
conda activate plink_env

python rss_correlation_of_borzoi_genetic_expression.py ${tissue_name} $borzoi_effect_size_file ${eqtl_sumstats_dir} ${protein_coding_genes_file} ${genotype_1000G_plink_stem} ${hm3_snp_list_file}