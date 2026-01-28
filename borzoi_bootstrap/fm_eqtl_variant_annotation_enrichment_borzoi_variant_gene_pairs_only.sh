#!/bin/bash
#SBATCH -t 0-4:50                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=20GB  




processed_fm_eqtl_output_file="${1}"
eqtl_sumstats_dir="${2}"
genotype_dir="${3}"
snp_anno_dir="${4}"
sig_thresh="${5}"
variant_annotation_enrichment_file="${6}"
remove_coding_boolean="${7}"
tissue_name="${8}"
borzoi_results_stem="${9}"
date
source ~/.bashrc
conda activate borzoi

python fm_eqtl_variant_annotation_enrichment_borzoi_variant_gene_pairs_only.py $processed_fm_eqtl_output_file $eqtl_sumstats_dir $genotype_dir $snp_anno_dir $sig_thresh $variant_annotation_enrichment_file ${remove_coding_boolean} ${tissue_name} $borzoi_results_stem
date