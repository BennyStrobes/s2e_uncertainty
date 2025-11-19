#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-7:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu                           # Partition to run in
#SBATCH --mem=7GB                         # Memory total in MiB (for all cores)

source ~/.bashrc
conda activate SAGEnet

fasta_file="${1}"
expr_file="${2}"
tss_data_file="${3}"
protein_coding_gene_list="${4}"
model_training_dir="${5}"
output_stem="${6}"

date

python train_rsage_on_gtex.py $fasta_file $expr_file $tss_data_file $protein_coding_gene_list $model_training_dir $output_stem

date