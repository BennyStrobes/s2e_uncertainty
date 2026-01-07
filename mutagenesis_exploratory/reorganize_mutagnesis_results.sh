#!/bin/bash
#SBATCH -t 0-5:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=4GB  





raw_mutagenesis_file="${1}"
reorganized_mutagenesis_file="${2}"
mutagenesis_variant_vcf="${3}"

source ~/.bashrc
conda activate borzoi

python reorganize_mutagenesis_results.py $raw_mutagenesis_file $reorganized_mutagenesis_file $mutagenesis_variant_vcf