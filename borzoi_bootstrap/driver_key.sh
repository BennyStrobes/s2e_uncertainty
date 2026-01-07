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

protein_coding_genes_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gene_tss.bed"

genotype_1000G_plink_stem="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/1000G.EUR.hg38."

hm3_snp_list_file="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/w_hm3.noMHC.snplist"

ldsc_snp_annotation_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

ldsc_hg19_cell_type_annotation_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3/cell_type_groups/"

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

# Organized borzoi gtex predictions
organized_borzoi_gtex_predictions=${output_root}"organized_borzoi_gtex_predictions/"

# bayes_posterior_effect_predictions
bayes_posterior_effect_predictions_dir=${output_root}"bayesian_borzoi_gtex_predictions/"

# Download borzoi data
borzoi_downloads_dir=${output_root}"borzoi_downloads/"

# Download borzoi data
expression_correlation_dir=${output_root}"expression_correlation/"

# Exploration of genome wide predictions dirctory
explore_genome_wide_pred_dir=${output_root}"genome_wide_pred_exploratory/"

# Variant annotation enrichments
variant_annotation_enrichment_dir=${output_root}"variant_annotation_enrichment/"

# Visualization dir
visualization_dir=${output_root}"visualize_borzoi/"

#################
# Run analysis
#################

##############################################################
# Get full target file (containing gtex tissue name too)
full_gtex_target_file="/lab-share/CHIP-Strober-e2/Public/Borzoi_data/w5/gtex_targets_full.txt"
gtex_tissue_names_file="/lab-share/CHIP-Strober-e2/Public/Borzoi_data/w5/gtex_tissues.txt"
gtex_eqtl_tissue_names_file="/lab-share/CHIP-Strober-e2/Public/Borzoi_data/w5/gtex_eqtl_tissues.txt"
if false; then
source ~/.bashrc
conda activate borzoi
python add_gtex_tissue_names_to_target_files.py $gtex_target_file $full_gtex_target_file $gtex_sample_attributes_file $gtex_tissue_names_file $gtex_eqtl_tissue_names_file $eqtl_sumstats_dir
fi

##############################################################
# Make TF-records 
# This is basically the Makefile in the borzoi tutorial (https://github.com/calico/borzoi/blob/main/tutorials/latest/make_data/Makefile)
# But in a shell script because i'm not fancy
# Has 75 folds built in
if false; then
sbatch make_tf_records.sh ${gtex_target_file} ${tfrecords_dir} ${baskerville_code_dir}
fi

##############################################################
# Prepare training/validation/test data splits
# Make file containing lines of submission instructions
submission_file=${model_training_dir}"submit_bootstrapped_gtex.txt"
if false; then
source ~/.bashrc
conda activate borzoi
python make_westminster_bootstrapped_train_folds.py -e borzoi --boot_start 1 --boot_end 100 -o ${model_training_dir}bootstrapped_models ${borzoi_micro_json_file} $submission_file ${tfrecords_dir}"hg38/"
fi


##############################################################
# Train models
submission_file=${model_training_dir}"submit_bootstrapped_gtex_first_10.txt"
# Submit training jobs
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi



# done
submission_file=${model_training_dir}"submit_bootstrapped_gtex_21_to_40.txt"
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi


# done
submission_file=${model_training_dir}"submit_bootstrapped_gtex_41_50.txt"
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi

# Running on code temp9
# Messed up: 72, 71, 70, 69, 67, 66, 54 ????
submission_file=${model_training_dirc}"submit_bootstrapped_gtex_51_75.txt"
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi


# Run on code_temp1
# messed up: 79
submission_file=${model_training_dir}"submit_bootstrapped_gtex_76_100.txt"
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi

# Running on code_temp2
if false; then
submission_file=${model_training_dir}"submit_bootstrapped_gtex_92_100.txt"
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi


# code_temp4
submission_file=${model_training_dir}"submit_bootstrapped_gtex_temper.txt"
if false; then
while IFS= read -r command; do
    sbatch train_borzoi.sh "$baskerville_code_dir$command"
done < "$submission_file"
fi






##############################################################
# Extract fine-mapped eqtls
pip_threshold="0.9"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"
fm_vcf_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl.vcf"
if false; then
sh parse_eqtl_data.sh $raw_fine_mapping_eqtl_results_file $pip_threshold $processed_fm_eqtl_output_file $fm_vcf_output_file
fi



##############################################################
# Get predicted Borzoi effects for fine-mapped eqtls
# Submit parallel jobs


# code_temp8
if false; then
for bs_iter in {1..100}; do
    output_file=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results.txt"
    sbatch fast_borzoi_sed.sh ${output_file} ${fm_vcf_output_file} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi

# Organize results across parallel runs
if false; then
sh organize_gtex_borzoi_sed_results_fm.sh ${borzoi_eqtl_effects_dir} $full_gtex_target_file $borzoi_eqtl_effects_dir $processed_fm_eqtl_output_file
fi


##############################################################
# Visualize fine-mapped eQTL analysis
if false; then
source ~/.bashrc
conda activate borzoi
python visualize_borzoi_estimated_effects.py $borzoi_eqtl_effects_dir"cross_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects_cross_tissue.txt" $visualization_dir
source ~/.bashrc
conda activate plink_env
Rscript visualize_borzoi_fine_mapped_effects.R $borzoi_eqtl_effects_dir"cross_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects_cross_tissue.txt" $visualization_dir
fi



##############################################################
# GO beyond just fine-mapped eqtls
# Generate list of variants in GTEx with MAF > 0.05 (in at least one tissue)
gtex_all_snps_vcf=${gtex_snp_dir}"gtex_snps.vcf"
if false; then
sbatch generate_vcf_for_all_gtex_snps.sh $eqtl_sumstats_dir $gtex_all_snps_vcf
fi





##############################################################
# Extract borzoi predicted effects for all gtex snps
part="0"
gtex_all_snps_vcf=${gtex_snp_dir}"gtex_snps.vcf.part.0"$part

# 70 hours of runtime here
if false; then
for bs_iter in {1..100}; do
    output_file=${borzoi_gtex_predictions}"bs"${bs_iter}"_part_"${part}"_borzoi_pred_eqtl_effects_borzoi_sed_results.txt"
    sbatch fast_borzoi_sed.sh ${output_file} ${gtex_all_snps_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi

if false; then
part="1"
gtex_all_snps_vcf=${gtex_snp_dir}"gtex_snps.vcf.part.0"$part
for bs_iter in {51..75}; do
    output_file=${borzoi_gtex_predictions}"bs"${bs_iter}"_part_"${part}"_borzoi_pred_eqtl_effects_borzoi_sed_results.txt"
    sbatch fast_borzoi_sed.sh ${output_file} ${gtex_all_snps_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi



# Organize results across parallel runs
if false; then
sbatch organize_gtex_borzoi_sed_results.sh ${borzoi_gtex_predictions} $full_gtex_target_file $organized_borzoi_gtex_predictions
fi


tissue_name="Muscle_Skeletal"
bayes_fit_output_root=${bayes_posterior_effect_predictions_dir}"gibbs_ash_model_"${tissue_name}"results"
if false; then
sh fit_bayesian_posterior_to_borzoi_effects.sh $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt" $bayes_fit_output_root
fi


tissue_name="Muscle_Skeletal"
expr_correlation_summary_file=${expression_correlation_dir}${tissue_name}"_expression_correlation_summary.txt"
if false; then
sh rss_correlation_of_borzoi_genetic_expression.sh $tissue_name $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt" $eqtl_sumstats_dir $protein_coding_genes_file $genotype_1000G_plink_stem $hm3_snp_list_file $expr_correlation_summary_file
fi



# Need to run for all tissues
if false; then
tail -n +2 $gtex_tissue_names_file | while read -r tissue_name; do
    sig_thresh="0.05"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_enrichments_"${tissue_name}"_"${sig_thresh}"_summary.txt"
    remove_coding="False"
    sbatch variant_annotation_enrichment.sh $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt" $genotype_1000G_plink_stem $ldsc_snp_annotation_dir $sig_thresh $variant_annotation_enrichment_file $remove_coding

    sig_thresh="0.05"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_enrichments_"${tissue_name}"_"${sig_thresh}"_remove_coding_summary.txt"
    remove_coding="True"
    sbatch variant_annotation_enrichment.sh $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt" $genotype_1000G_plink_stem $ldsc_snp_annotation_dir $sig_thresh $variant_annotation_enrichment_file $remove_coding
done
fi


# Minimal signal expected
if false; then
tail -n +2 $gtex_tissue_names_file | while read -r tissue_name; do
    echo "Processing ${tissue_name}"

    sig_thresh="0.05"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_cell_type_enrichments_"${tissue_name}"_"${sig_thresh}"_summary.txt"
    sbatch variant_cell_type_annotation_enrichment.sh $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt" $genotype_1000G_plink_stem $ldsc_hg19_cell_type_annotation_dir $sig_thresh $variant_annotation_enrichment_file

done
fi







if false; then
source ~/.bashrc
conda activate borzoi
tissue_name="Muscle_Skeletal"
python explore_genome_wide_predictions.py $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt" $explore_genome_wide_pred_dir
fi


















# OLD

# Organize results across parallel jobs
if false; then
sh organize_borzoi_predicted_eqtl_effects.sh ${borzoi_eqtl_effects_dir} $pip_threshold $processed_fm_eqtl_output_file $gtex_sample_attributes_file ${model_training_dir}"bootstrapped_models/"
fi


