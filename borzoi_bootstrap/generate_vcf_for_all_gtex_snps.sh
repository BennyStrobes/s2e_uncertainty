#!/bin/bash
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB 




gtex_sumstat_dir="${1}"
output_vcf_file="${2}"

source ~/.bashrc
conda activate borzoi


python generate_vcf_for_all_gtex_snps.py $gtex_sumstat_dir $output_vcf_file