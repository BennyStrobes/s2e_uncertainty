#!/bin/bash
#SBATCH -t 0-6:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=30GB 



borzoi_effect_file="${1}"
eqtl_effects_file="${2}"
borzoi_annotation_file="${3}"
genotype_stem="${4}"
genotype_sample_mapping_file="${5}"
ld_corr_output_stem="${6}"



source ~/.bashrc
conda activate plink_env







if false; then
python run_ld_corr.py $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $ld_corr_output_stem
fi

python run_ld_corr_HE.py $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $ld_corr_output_stem
