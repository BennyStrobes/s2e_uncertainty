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



python "borzoi_sed.py" -o ${output_dir}"rc_only" --rc --stats logSED,logD2,refLogSed,altLogSed -t ${borzoi_training_dir}"data0/targets.txt" ${borzoi_training_dir}"params.json" ${borzoi_training_dir}"train/model_best.h5" $vcf_input_file

if false; then


python "borzoi_sed.py" -o ${output_dir}"rc_only_5_shifts" --rc --shifts "-2,-1,0,1,2" --stats logSED,logD2,refLogSed,altLogSed -t ${borzoi_training_dir}"data0/targets.txt" ${borzoi_training_dir}"params.json" ${borzoi_training_dir}"train/model_best.h5" $vcf_input_file


python "borzoi_sed.py" -o ${output_dir}"rc_only_9_shifts" --rc --shifts "-4,-3,-2,-1,0,1,2,3,4" --stats logSED,logD2,refLogSed,altLogSed -t ${borzoi_training_dir}"data0/targets.txt" ${borzoi_training_dir}"params.json" ${borzoi_training_dir}"train/model_best.h5" $vcf_input_file



python "borzoi_sed.py" -o ${output_dir}"rc_only_9_big_shifts" --rc --shifts "-15,-10,-5,0,5,10,15" --stats logSED,logD2,refLogSed,altLogSed -t ${borzoi_training_dir}"data0/targets.txt" ${borzoi_training_dir}"params.json" ${borzoi_training_dir}"train/model_best.h5" $vcf_input_file
fi