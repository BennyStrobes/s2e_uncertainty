#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-3:20                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu-pe                           # Partition to run in
#SBATCH --mem=10GB  



source ~/.bashrc
conda activate borzoi
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

output_dir="${1}"
vcf_input_file="${2}"
borzoi_training_dir="${3}"
processed_fm_eqtl_output_file="${4}"
full_gtex_target_file="${5}"
gene_tss_file="${6}"

echo $output_dir

date
python "fast_borzoi_sed_per_bin.py" -o ${output_dir} --rc --stats logSED,refLog,altLog -t ${borzoi_training_dir}"data0/targets.txt" -z $processed_fm_eqtl_output_file -a ${full_gtex_target_file} -x ${gene_tss_file} ${borzoi_training_dir}"params.json" ${borzoi_training_dir}"train/model_best.h5" $vcf_input_file
date