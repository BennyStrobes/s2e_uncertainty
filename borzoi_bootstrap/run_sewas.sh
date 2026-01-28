#!/bin/bash
#SBATCH -t 0-4:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=10GB  



trait_sumstat_file="${1}"
borzoi_results_file="${2}"
genotype_1000G_plink_stem="${3}"
sewas_output_file="${4}"
protein_coding_genes_file="${5}"

source ~/.bashrc
conda activate plink_env


python SeWAS.py \
	--sumstat-file $trait_sumstat_file \
	--borzoi-result-file $borzoi_results_file \
	--plink-genotype-stem $genotype_1000G_plink_stem \
	--output-file ${sewas_output_file} \
	--gene-list-file ${protein_coding_genes_file}