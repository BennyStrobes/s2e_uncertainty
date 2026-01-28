#!/bin/bash
#SBATCH -t 0-15:50                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=20GB  




source ~/.bashrc
conda activate borzoi


borzoi_results_stem="${1}"
processed_fm_eqtl_output_file="${2}"
eqtl_sumstats_dir="${3}"
gtex_eqtl_tissue_names_file="${4}"
tissue_name="${5}"
sig_thresh="${6}"
variant_annotation_enrichment_file="${7}"


echo $tissue_name

python borzoi_variant_enrichment_in_fm_eqtls.py $borzoi_results_stem $processed_fm_eqtl_output_file $eqtl_sumstats_dir $gtex_eqtl_tissue_names_file $tissue_name $sig_thresh $variant_annotation_enrichment_file
