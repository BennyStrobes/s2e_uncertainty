#!/bin/bash
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 



borzoi_standardized_tracks_file="${1}"
eqtl_effects_file="${2}"
genotype_stem="${3}"
ld_corr_compete_output_stem="${4}"
model="${5}"
borzoi_thresh="${6}"


source ~/.bashrc
conda activate plink_env

echo $ld_corr_compete_output_stem


python run_ld_corr_compete.py $borzoi_standardized_tracks_file $eqtl_effects_file $genotype_stem $ld_corr_compete_output_stem $model $borzoi_thresh