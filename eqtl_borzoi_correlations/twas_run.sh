#!/bin/bash
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 





borzoi_effect_file="${1}"
borzoi_annotation_file="${2}"
genotype_stem="${3}"
genotype_sample_mapping_file="${4}"
expr_file="${5}"
distribution="${6}"
borzoi_based_prior_output_stem="${7}"
twas_output_file="${8}"
standardize_geno="${9}"
trait_ss_file="${10}"

source ~/.bashrc
conda activate plink_env

date

echo $twas_output_file

python twas_run.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $twas_output_file $standardize_geno $trait_ss_file
date