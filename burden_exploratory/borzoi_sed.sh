#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-13:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu                           # Partition to run in
#SBATCH --mem=10GB  



source ~/.bashrc
conda activate borzoi
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

output_dir="${1}"
vcf_input_file="${2}"
borzoi_training_dir="${3}"
model_num="${4}"
split_number="${5}"


echo $model_num"_"${split_number}

python "borzoi_sed.py" -o ${output_dir} --rc --stats logSED -t ${borzoi_training_dir}"targets_gtex.txt" ${borzoi_training_dir}"params_pred.json" ${borzoi_training_dir}"model0_best_f3c"${model_num}".h5" $vcf_input_file
