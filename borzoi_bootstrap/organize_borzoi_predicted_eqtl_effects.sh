#!/bin/bash
#SBATCH -t 0-2:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=10GB  


borzoi_eqtl_effects_dir="${1}"
pip_threshold="${2}"
processed_fm_eqtl_output_file="${3}"
gtex_sample_attributes_file="${4}"
model_training_dir="${5}"


source ~/.bashrc
conda activate borzoi

for bs_iter in {1..50}; do


	echo "Bootstrap: "${bs_iter}
	borzoi_eqtl_output_dir=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results"
	borzoi_eqtl_output_file=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_cross_tissue.txt"

	python extract_borzoi_effects_cross_tissues.py $processed_fm_eqtl_output_file ${borzoi_eqtl_output_dir}"/sed.h5" $borzoi_eqtl_output_file $gtex_sample_attributes_file ${model_training_dir}"bs"${bs_iter}"/data0/targets.txt"
done




organized_bs_eqtl_output=${borzoi_eqtl_effects_dir}"cross_bootstrap_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_cross_tissue.txt"
python organize_borzoi_predicted_effects_across_bootstraps.py "1" "50" ${borzoi_eqtl_effects_dir} ${organized_bs_eqtl_output}
