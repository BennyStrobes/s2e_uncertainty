#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=40GB  





borzoi_results_file_stem="${1}"
processed_fm_eqtl_output_file="${2}"
eqtl_sumstats_dir="${3}"
gene_tss_file="${4}"
dist_to_tss_summary_file="${5}"
tissue_name="${6}"


source ~/.bashrc
conda activate borzoi



python dist_to_tss_enrichment.py $borzoi_results_file_stem $processed_fm_eqtl_output_file $eqtl_sumstats_dir $gene_tss_file $dist_to_tss_summary_file $tissue_name
