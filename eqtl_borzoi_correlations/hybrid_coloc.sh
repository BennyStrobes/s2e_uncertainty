#!/bin/bash
#SBATCH -t 0-14:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=30GB 





borzoi_effect_file="${1}"
borzoi_annotation_file="${2}"
genotype_stem="${3}"
genotype_sample_mapping_file="${4}"
expr_file="${5}"
coloc_output_file="${6}"
standardize_geno="${7}"
training_sample_size="${8}"
trait_ss_file="${9}"

source ~/.bashrc
conda activate plink_env

echo $coloc_output_file

python hybrid_coloc.py $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $coloc_output_file $standardize_geno $training_sample_size $trait_ss_file