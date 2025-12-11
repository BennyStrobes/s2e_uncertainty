#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-10:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)

#################
# Input data
#################
gtex_target_file="/lab-share/CHIP-Strober-e2/Public/Borzoi_data/w5/gtex_targets.txt"

baskerville_code_dir="/home/ch271704/tools/baskerville/src/baskerville/scripts/"

borzoi_micro_json_file="/home/ch271704/tools/borzoi/tutorials/latest/train_model/params_gtex_micro.json"
borzoi_micro_json_file="/home/ch271704/tools/borzoi/tutorials/latest/train_model/params_gtex_micro_100_iter.json"

raw_fine_mapping_eqtl_results_file="/lab-share/CHIP-Strober-e2/Public/GTEx/fine_mapping/GTEx_49tissues_release1.tsv"

gtex_sample_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/gtex_sample_attributes/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"

eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/eqtl_sumstats/"

borzoi_code_dir="/home/ch271704/tools/borzoi/src/scripts/"
#################
# Output directories
#################
# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_tissue_bootstrap/"

# Directory containing processed tfrecords
tfrecords_dir=${output_root}"tfrecords/"

# Directory containing model training results
model_training_dir=${output_root}"model_train/"

# Directory containing processed Fine-mapped eqtl results
processed_fm_eqtl_data_dir=${output_root}"processed_fm_eqtl_results/"

# Directory contain borzoi predicted fine-mapped eqtl results
borzoi_eqtl_effects_dir=${output_root}"borzoi_pred_eqtl_effects/"

# GTEx snp dir
gtex_snp_dir=${output_root}"GTEx_snps/"

# Borzoi gtex predictions
borzoi_gtex_predictions=${output_root}"borzoi_gtex_predicted_effects/"

# Download borzoi data
borzoi_downloads_dir=${output_root}"borzoi_downloads/"

# Visualization dir
visualization_dir=${output_root}"visualize_borzoi/"

#################
# Run analysis
#################

# Make TF-records 
# This is basically the Makefile in the borzoi tutorial (https://github.com/calico/borzoi/blob/main/tutorials/latest/make_data/Makefile)
# But in a shell script because i'm not fancy
# Has 75 folds built in
if false; then
sbatch make_tf_records.sh ${gtex_target_file} ${tfrecords_dir} ${baskerville_code_dir}
fi

# Prepare training/validation/test data splits
# Make file containing lines of submission instructions
submission_file=${model_training_dir}"submit_bootstrapped_gtex.txt"
if false; then
source ~/.bashrc
conda activate borzoi
python make_westminster_bootstrapped_train_folds.py -e borzoi --boot_start 1 --boot_end 100 -o ${model_training_dir}bootstrapped_models ${borzoi_micro_json_file} $submission_file ${tfrecords_dir}"hg38/"
fi

submission_file=${model_training_dir}"submit_bootstrapped_gtex_first_10.txt"
# Submit training jobs
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi



# Running on code_temp6
submission_file=${model_training_dir}"submit_bootstrapped_gtex_21_to_40.txt"
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi


# Run on code_temp9
submission_file=${model_training_dir}"submit_bootstrapped_gtex_41_50.txt"
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi



pip_threshold="0.9"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"
fm_vcf_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl.vcf"
if false; then
sh parse_eqtl_data.sh $raw_fine_mapping_eqtl_results_file $pip_threshold $processed_fm_eqtl_output_file $fm_vcf_output_file
fi


if false; then
for bs_iter in {1..20}; do
    borzoi_eqtl_output_dir=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results"
    sbatch borzoi_sed.sh ${borzoi_eqtl_output_dir} ${fm_vcf_output_file} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi

if false; then
for bs_iter in {41..50}; do
    borzoi_eqtl_output_dir=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results"
    sbatch borzoi_sed.sh ${borzoi_eqtl_output_dir} ${fm_vcf_output_file} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi

if false; then
sh organize_borzoi_predicted_eqtl_effects.sh ${borzoi_eqtl_effects_dir} $pip_threshold $processed_fm_eqtl_output_file $gtex_sample_attributes_file ${model_training_dir}"bootstrapped_models/"
fi

if false; then
source ~/.bashrc
conda activate borzoi
python visualize_borzoi_estimated_effects.py $borzoi_eqtl_effects_dir"cross_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects_cross_tissue.txt" $visualization_dir
fi




# GO beyond just fine-mapped eqtls
gtex_all_snps_vcf=${gtex_snp_dir}"gtex_snps.vcf"
if false; then
sh generate_vcf_for_all_gtex_snps.sh $eqtl_sumstats_dir $gtex_all_snps_vcf
fi



gtex_all_snps_vcf=${gtex_snp_dir}"gtex_snps_short.vcf"

bs_iter="1"
output_file=${borzoi_gtex_predictions}"bs"${bs_iter}"_borzoi_pred_eqtl_effects_borzoi_sed_results.txt"
if false; then
sbatch fast_borzoi_sed.sh ${output_file} ${fm_vcf_output_file} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
fi



















# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/old_gtex_tissue_bootstrap/"

# Directory containing processed tfrecords
tfrecords_dir=${output_root}"tfrecords/"

# Directory containing model training results
model_training_dir=${output_root}"model_train/"

# Directory containing processed Fine-mapped eqtl results
processed_fm_eqtl_data_dir=${output_root}"processed_fm_eqtl_results/"

# Directory contain borzoi predicted fine-mapped eqtl results
borzoi_eqtl_effects_dir=${output_root}"borzoi_pred_eqtl_effects/"

# Download borzoi data
borzoi_downloads_dir=${output_root}"borzoi_downloads/"

# Visualization dir
visualization_dir=${output_root}"visualize_borzoi/"


if false; then
source ~/.bashrc
conda activate borzoi

python visualize_borzoi_estimated_effects.py $borzoi_eqtl_effects_dir"cross_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects_cross_tissue.txt" $visualization_dir
fi













#################
# OLD
#################






# Download borzoi data
if false; then
sh download_borzoi_data.sh $borzoi_downloads_dir
fi


pip_threshold="0.9"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"
if false; then
sh parse_eqtl_data.sh $raw_fine_mapping_eqtl_results_file $pip_threshold $processed_fm_eqtl_output_file
fi

fm_vcf_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl.vcf"
if false; then
source ~/.bashrc
conda activate borzoi
python make_vcf_out_of_eqtls.py $processed_fm_eqtl_output_file $fm_vcf_output_file
fi


borzoi_eqtl_output_dir=${borzoi_eqtl_effects_dir}"PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results"
if false; then
source ~/.bashrc
conda activate borzoi
python "/home/ch271704/tools/borzoi/src/scripts/borzoi_sed.py" -o ${borzoi_eqtl_output_dir} --rc --stats logSED,logD2 -t ${model_training_dir}"micro_models/f0c0/data0/targets.txt" ${model_training_dir}"micro_models/f0c0/params.json" ${model_training_dir}"micro_models/f0c0/train/model_best.h5" $fm_vcf_output_file
fi

if false; then
borzoi_eqtl_output_file=${borzoi_eqtl_effects_dir}"PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_cross_tissue.txt"
python extract_borzoi_effects_cross_tissues.py $processed_fm_eqtl_output_file ${borzoi_eqtl_output_dir}"/sed.h5" $borzoi_eqtl_output_file $gtex_sample_attributes_file ${model_training_dir}"micro_models/f0c0/data0/targets.txt"
fi















##############
# OLD
##############





#################
# Evaluate training using fine-mapped eqtl effect sizes
# Parse eQTL data
tissue_name="Whole_Blood"
pip_threshold="0.9"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}$tissue_name"_PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"
if false; then
sh parse_eqtl_data.sh $raw_fine_mapping_eqtl_results_file $tissue_name $pip_threshold $processed_fm_eqtl_output_file
fi



# Predict fine-mapped eqtl effect sizes with borzoi
if false; then
    sim_iter="0"
    total_sims="1"
    borzoi_eqtl_output_file=${borzoi_eqtl_effects_dir}${tissue_name}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_"$sim_iter"_"${total_sims}".txt"
    sbatch extract_borzoi_effects.sh $processed_fm_eqtl_output_file $borzoi_downloads_dir $model_training_dir $borzoi_eqtl_output_file $sim_iter $total_sims

fi





fm_vcf_output_file=${processed_fm_eqtl_data_dir}$tissue_name"_PIP_"${pip_threshold}"_fine_mapped_eqtl.vcf"
if false; then
python make_vcv_out_of_eqtls.py $processed_fm_eqtl_output_file $fm_vcf_output_file
fi

borzoi_eqtl_output_file2=${borzoi_eqtl_effects_dir}${tissue_name}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_"$sim_iter"_"${total_sims}"_v2.txt"
if false; then
python "/home/ch271704/tools/borzoi/src/scripts/borzoi_sed.py" -o ${borzoi_eqtl_output_file2} --rc --stats logSED,logD2 -t ${model_training_dir}"micro_models/f0c0/data0/targets.txt" ${model_training_dir}"micro_models/f0c0/params.json" ${model_training_dir}"micro_models/f0c0/train/model_best.h5" $fm_vcf_output_file
fi

borzoi_eqtl_output_file3=${borzoi_eqtl_effects_dir}${tissue_name}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_"$sim_iter"_"${total_sims}"_v3.txt"
if false; then
python extract_borzoi_effects2.py $processed_fm_eqtl_output_file ${borzoi_eqtl_output_file2}"/sed.h5" $borzoi_eqtl_output_file3
fi



