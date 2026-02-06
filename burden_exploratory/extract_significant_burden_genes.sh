#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)



trait_names_file="${1}"
burden_test_data_dir="${2}"
sig_burden_genes_dir="${3}"
gene_annotation_summary_file="${4}"


source ~/.bash_profile

python extract_significant_burden_genes.py $trait_names_file $burden_test_data_dir $sig_burden_genes_dir $gene_annotation_summary_file
