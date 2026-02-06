#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-25:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)







burden_genes_summary_file="${1}"
genotype_dir="${2}"
variant_gene_pairs_dir="${3}"

source ~/.bash_profile


python extract_variant_gene_pairs_to_test.py $burden_genes_summary_file $genotype_dir $variant_gene_pairs_dir