#!/bin/bash
#SBATCH -t 0-5:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=20GB  




borzoi_effect_est_file="${1}"
genotype_dir="${2}"
snp_anno_dir="${3}"
sig_thresh="${4}"
variant_annotation_enrichment_file="${5}"
remove_coding_boolean="${6}"

date
source ~/.bashrc
conda activate borzoi

python variant_annotation_enrichment.py $borzoi_effect_est_file $genotype_dir $snp_anno_dir $sig_thresh $variant_annotation_enrichment_file ${remove_coding_boolean}
date