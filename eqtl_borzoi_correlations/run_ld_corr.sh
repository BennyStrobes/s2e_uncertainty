#!/bin/bash
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 



borzoi_effect_file="${1}"
eqtl_effects_file="${2}"
borzoi_annotation_file="${3}"
genotype_stem="${4}"
ld_corr_output_stem="${5}"
chi_sq_thresh="${6}"


source ~/.bashrc
conda activate plink_env








python run_ld_corr.py $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem $chi_sq_thresh
