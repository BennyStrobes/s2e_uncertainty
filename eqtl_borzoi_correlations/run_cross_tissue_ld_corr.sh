#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=40GB 



t1_borzoi_effect_file="${1}"
t2_borzoi_effect_file="${2}"
t1_eqtl_effects_file="${3}"
t2_eqtl_effects_file="${4}"
t1_borzoi_annotation_file="${5}"
t2_borzoi_annotation_file="${6}"
genotype_stem="${7}"
t1_genotype_sample_mapping_file="${8}"
t2_genotype_sample_mapping_file="${9}"
cross_tissue_ld_corr_output_stem="${10}"




source ~/.bashrc
conda activate plink_env

echo $cross_tissue_ld_corr_output_stem



python run_cross_tissue_ld_corr.py $t1_borzoi_effect_file $t2_borzoi_effect_file $t1_eqtl_effects_file $t2_eqtl_effects_file $t1_borzoi_annotation_file $t2_borzoi_annotation_file $genotype_stem $t1_genotype_sample_mapping_file $t2_genotype_sample_mapping_file $cross_tissue_ld_corr_output_stem
