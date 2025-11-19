#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-1:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu                           # Partition to run in
#SBATCH --mem=5GB                         # Memory total in MiB (for all cores)

source ~/.bashrc
conda activate SAGEnet






processed_fm_eqtl_output_file="${1}"
hg38_fasta_file="${2}"
tss_data_file="${3}"
protein_coding_gene_list="${4}"
model_file="${5}"
output_file="${6}"
n_bootstraps="${7}"


python get_bs_delta_scores_for_fm_eqtl_variants.py $processed_fm_eqtl_output_file $hg38_fasta_file $tss_data_file $protein_coding_gene_list $model_file $output_file ${n_bootstraps}
