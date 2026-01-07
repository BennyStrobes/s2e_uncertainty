#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-0:20                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu-pe                           # Partition to run in
#SBATCH --mem=10GB  



source ~/.bashrc
conda activate borzoi
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

output_dir="${1}"
vcf_input_file="${2}"
borzoi_training_dir="${3}"

echo $output_dir

date
python "fast_borzoi_sed.py" -o ${output_dir} --rc --stats logSED,refLog,altLog -t ${borzoi_training_dir}"data0/targets.txt" ${borzoi_training_dir}"params.json" ${borzoi_training_dir}"train/model_best.h5" $vcf_input_file
date