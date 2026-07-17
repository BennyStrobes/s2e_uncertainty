#!/bin/bash
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



borzoi_effect_file="${1}"
annotation_name_file="${2}"
borzoi_annotation_file="${3}"
eqtl_sumstats_file="${4}"
tissue_name="${5}"
baselineLD_anno_dir="${6}"

source ~/.bashrc
conda activate plink_env


echo $borzoi_annotation_file

python annotate_variant_gene_pairs.py $borzoi_effect_file $annotation_name_file $borzoi_annotation_file $eqtl_sumstats_file $tissue_name $baselineLD_anno_dir
