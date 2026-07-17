#!/bin/bash
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=40GB 



borzoi_effect_file="${1}"
eqtl_effects_file="${2}"
borzoi_annotation_file="${3}"
genotype_stem="${4}"
genotype_sample_mapping_file="${5}"
bootstrapped_cross_tissue_gene_sets_file_stem="${6}"
ld_corr_output_stem="${7}"


source ~/.bashrc
conda activate plink_env

date

python run_ld_corr.py \
	--est-borzoi-effect-size-file $borzoi_effect_file \
	--est-eqtl-effect-size-file $eqtl_effects_file \
	--sim-variant-gene-annotation-file $borzoi_annotation_file \
	--genotype-plink-filestem $genotype_stem \
	--genotype-sample-mapping-file $genotype_sample_mapping_file \
	--ld-corr-output-stem $ld_corr_output_stem"_global_bs" \
	--bootstrapped-gene-set-filestem ${bootstrapped_cross_tissue_gene_sets_file_stem} \
    --weighted "True"


date

if false; then
# Without cross-tissue shared bootstraps
python run_ld_corr.py \
	--est-borzoi-effect-size-file $borzoi_effect_file \
	--est-eqtl-effect-size-file $eqtl_effects_file \
	--sim-variant-gene-annotation-file $borzoi_annotation_file \
	--genotype-plink-filestem $genotype_stem \
	--genotype-sample-mapping-file $genotype_sample_mapping_file \
	--ld-corr-output-stem $ld_corr_output_stem"_local_bs" \
    --weighted "True"
fi





###############################
# Old. not currently used
###############################

if false; then
python run_ld_corr_HE.py $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $ld_corr_output_stem
fi