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

raw_fine_mapping_eqtl_results_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/fine_mapping/v10/"

gtex_sample_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/gtex_sample_attributes/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"

eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/eqtl_sumstats/"

borzoi_code_dir="/home/ch271704/tools/borzoi/src/scripts/"

protein_coding_genes_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gene_tss.bed"

genotype_1000G_plink_stem="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/1000G.EUR.hg38."

hm3_snp_list_file="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/w_hm3.noMHC.snplist"

ldsc_snp_annotation_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

ldsc_weights_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

ldsc_hg19_cell_type_annotation_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3/cell_type_groups/"

ldsc_code_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/ldsc/"

gwas_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/sumstats_formatted_2024/sumstats/"

non_redundant_gwas_traits_file="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"

non_redundant_gwas_traits_file="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/sumstats_formatted_2024/non_redundent_traits_EUR_2024.txt"

gwas_sldsc_results_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/sumstats_formatted_2024/sldsc_h2/"

ukbb_all_snp_sumstat_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/UKBB_all_snps_sumstats/data/"

gene_tss_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/genecode.v39.GRCh38.bed"

gtex_v10_gene_gtf_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.GRCh38.genes.gtf"

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

sldsc_processed_anno_dir=${output_root}"sldsc_anno/"

sldsc_h2_results_dir=${output_root}"sldsc/"

# Visualization dir
visualization_fm_eqtl_dir=${output_root}"visualize_borzoi_fm_eqtls/"

# Visualize variant annotation enrichment
visualize_variant_anno_enrich_dir=${output_root}"visualize_variant_anno_enrichment/"

# Directory containing expression correlation visualizations
visualize_expression_correlation_dir=${output_root}"visualize_expression_correlation/"

# Directory containing SeWAS results
sewas_results_dir=${output_root}"SeWAS_results/"

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
sh parse_eqtl_data.sh $raw_fine_mapping_eqtl_results_dir $pip_threshold $processed_fm_eqtl_output_file $fm_vcf_output_file
fi

if false; then
pip_threshold="0.75"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"
fm_vcf_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl.vcf"
sh parse_eqtl_data.sh $raw_fine_mapping_eqtl_results_file $pip_threshold $processed_fm_eqtl_output_file $fm_vcf_output_file

pip_threshold="0.5"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"
fm_vcf_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_threshold}"_fine_mapped_eqtl.vcf"
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


if false; then
for bs_iter in {51..100}; do
    output_file=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_gene_span_borzoi_sed_results.txt"
    sbatch fast_borzoi_sed_gene_span.sh ${output_file} ${fm_vcf_output_file} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi



if false; then
for bs_iter in {31..60}; do
    output_file=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_window_span_borzoi_sed_results.txt"
    sbatch fast_borzoi_sed_window_span.sh ${output_file} ${fm_vcf_output_file} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi




tmp_vcf="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_tissue_bootstrap/processed_fm_eqtl_results/PIP_0.9_fine_mapped_eqtl_debug.vcf"
orig_gtf="/home/ch271704/tools/borzoi/examples/hg38/genes/gencode41/gencode41_basic_nort.gtf"
if false; then
for bs_iter in {1..100}; do
    output_file=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results_debug_fast_borzoi_pred.txt"
    sbatch fast_borzoi_sed_custom_gtf.sh ${output_file} ${tmp_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/" $orig_gtf
done
fi

if false; then
for bs_iter in {1..100}; do
    output_dir=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results_debug_"
    sbatch borzoi_sed_debugger.sh ${output_dir} ${tmp_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi

if false; then
for bs_iter in {1..30}; do
    output_dir=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results_debug_"
    sh borzoi_sed_debugger.sh ${output_dir} ${tmp_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi
if false; then
for bs_iter in {1..20}; do
    output_dir=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results_debug_"
    sh borzoi_sed_debugger.sh ${output_dir} ${tmp_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi


if false; then
bs_iter="2"
    output_dir=${borzoi_eqtl_effects_dir}"bs"${bs_iter}"_PIP_"${pip_threshold}"_borzoi_pred_eqtl_effects_borzoi_sed_results_debug_"
    sh borzoi_sed_debugger.sh ${output_dir} ${tmp_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
fi

# Organize results across parallel runs
if false; then
sh organize_gtex_borzoi_sed_results_fm.sh ${borzoi_eqtl_effects_dir} $full_gtex_target_file $borzoi_eqtl_effects_dir $processed_fm_eqtl_output_file $gene_tss_file
fi


##############################################################
# Visualize fine-mapped eQTL analysis
if false; then
source ~/.bashrc
conda activate borzoi
python visualize_borzoi_estimated_effects.py $borzoi_eqtl_effects_dir"cross_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects_cross_tissue.txt" $visualization_fm_eqtl_dir
fi

if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_borzoi_fine_mapped_effects.R $borzoi_eqtl_effects_dir"cross_bootstrap_PIP_0.9_gene_span_borzoi_pred_eqtl_effects_cross_tissue.txt" $visualization_fm_eqtl_dir
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
for bs_iter in {76..100}; do
    output_file=${borzoi_gtex_predictions}"bs"${bs_iter}"_part_"${part}"_borzoi_pred_eqtl_effects_borzoi_sed_results.txt"
    sbatch fast_borzoi_sed.sh ${output_file} ${gtex_all_snps_vcf} ${model_training_dir}"bootstrapped_models/bs"${bs_iter}"/"
done
fi


###################
# Organize results across parallel runs
part="1"
if false; then
sbatch organize_gtex_borzoi_sed_results.sh ${borzoi_gtex_predictions} $full_gtex_target_file $organized_borzoi_gtex_predictions $part
fi


###################
# Fit bayesian posterior to results (doesn't seem to be working)
tissue_name="Muscle_Skeletal"
bayes_fit_output_root=${bayes_posterior_effect_predictions_dir}"gibbs_ash_model_"${tissue_name}"results"
if false; then
sh fit_bayesian_posterior_to_borzoi_effects.sh $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt" $bayes_fit_output_root
fi

###################
# Compute correlations with predicted expression
tissue_name="Muscle_Skeletal"
expr_correlation_summary_file=${expression_correlation_dir}${tissue_name}"_expression_correlation_summary.txt"
if false; then
sh rss_correlation_of_borzoi_genetic_expression.sh $tissue_name $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_" $eqtl_sumstats_dir $protein_coding_genes_file $genotype_1000G_plink_stem $hm3_snp_list_file $expr_correlation_summary_file
fi



if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_expression_correlations.R $tissue_name $expr_correlation_summary_file $visualize_expression_correlation_dir
fi



#####################
# Run SeWAS
tissue_name="Liver"
trait_name="biochemistry_Cholesterol"
trait_sumstat_file=${ukbb_all_snp_sumstat_dir}${trait_name}"_hg38_liftover_sumstats.bgen.stats.gz"
borzoi_results_file=$organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt"

sewas_output_file=${sewas_results_dir}"SeWAS_"${trait_name}"_"${tissue_name}"_all_protein_coding_genes.txt"
if false; then
sh run_sewas.sh ${trait_sumstat_file} ${borzoi_results_file} $genotype_1000G_plink_stem ${sewas_output_file} $protein_coding_genes_file
fi







###################
# Assess which annotations are enriched in significant borzoi effects
# Need to run for all tissues
if false; then
tail -n +2 $gtex_tissue_names_file | while read -r tissue_name; do

    sig_thresh="0.05"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_enrichments_"${tissue_name}"_"${sig_thresh}"_summary.txt"
    remove_coding="False"
    sbatch variant_annotation_enrichment.sh $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_" $genotype_1000G_plink_stem $ldsc_snp_annotation_dir $sig_thresh $variant_annotation_enrichment_file $remove_coding $eqtl_sumstats_dir $tissue_name

    sig_thresh="0.05"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_enrichments_"${tissue_name}"_"${sig_thresh}"_remove_coding_summary.txt"
    remove_coding="True"
    sbatch variant_annotation_enrichment.sh $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_" $genotype_1000G_plink_stem $ldsc_snp_annotation_dir $sig_thresh $variant_annotation_enrichment_file $remove_coding $eqtl_sumstats_dir $tissue_name
done
fi







tissue_name="Whole_Blood"
pip_thresh="0.5"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_thresh}"_fine_mapped_eqtl_results.txt"
dist_to_tss_summary_file=$variant_annotation_enrichment_dir"dist_to_tss_summary_"${tissue_name}".txt"
if false; then
sbatch dist_to_tss_enrichment.sh $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_" $processed_fm_eqtl_output_file $eqtl_sumstats_dir $gene_tss_file $dist_to_tss_summary_file $tissue_name
fi



# Run variant enrichment for fine-mapped eQTLs
if false; then
tail -n +2 $gtex_eqtl_tissue_names_file | while read -r tissue_name; do
    pip_thresh="0.9"
    processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_thresh}"_fine_mapped_eqtl_results.txt"
    eqtl_variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"fm_eqtl_variant_enrichments_"${tissue_name}"_"${pip_thresh}"_summary.txt"
    remove_coding="False"
    sbatch fm_eqtl_variant_annotation_enrichment.sh $processed_fm_eqtl_output_file $eqtl_sumstats_dir $genotype_1000G_plink_stem $ldsc_snp_annotation_dir $pip_thresh $eqtl_variant_annotation_enrichment_file $remove_coding $tissue_name
    

    pip_thresh="0.9"
    processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_thresh}"_fine_mapped_eqtl_results.txt"
    eqtl_variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"fm_eqtl_borzoi_vg_pairs_only_variant_enrichments_"${tissue_name}"_"${pip_thresh}"_summary.txt"
    remove_coding="False"
    sbatch fm_eqtl_variant_annotation_enrichment_borzoi_variant_gene_pairs_only.sh $processed_fm_eqtl_output_file $eqtl_sumstats_dir $genotype_1000G_plink_stem $ldsc_snp_annotation_dir $pip_thresh $eqtl_variant_annotation_enrichment_file $remove_coding $tissue_name $organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_"
done
fi


if false; then
tail -n +2 $gtex_eqtl_tissue_names_file | while read -r tissue_name; do
    pip_thresh="0.9"
    sig_thresh="0.05"
    processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_thresh}"_fine_mapped_eqtl_results.txt"
    borzoi_results_stem=$organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_enrichment_in_fine_mapped_eqtls_"${tissue_name}"_"${sig_thresh}"_"${pip_thresh}"_summary.txt"
    sbatch borzoi_variant_enrichment_in_fm_eqtls.sh $borzoi_results_stem $processed_fm_eqtl_output_file $eqtl_sumstats_dir $gtex_eqtl_tissue_names_file $tissue_name $sig_thresh $variant_annotation_enrichment_file


    pip_thresh="0.9"
    sig_thresh="0.15"
    processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_thresh}"_fine_mapped_eqtl_results.txt"
    borzoi_results_stem=$organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_enrichment_in_fine_mapped_eqtls_"${tissue_name}"_"${sig_thresh}"_"${pip_thresh}"_summary.txt"
    sbatch borzoi_variant_enrichment_in_fm_eqtls.sh $borzoi_results_stem $processed_fm_eqtl_output_file $eqtl_sumstats_dir $gtex_eqtl_tissue_names_file $tissue_name $sig_thresh $variant_annotation_enrichment_file

    pip_thresh="0.5"
    sig_thresh="0.2"
    processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_thresh}"_fine_mapped_eqtl_results.txt"
    borzoi_results_stem=$organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_enrichment_in_fine_mapped_eqtls_"${tissue_name}"_"${sig_thresh}"_"${pip_thresh}"_summary.txt"
    sbatch borzoi_variant_enrichment_in_fm_eqtls.sh $borzoi_results_stem $processed_fm_eqtl_output_file $eqtl_sumstats_dir $gtex_eqtl_tissue_names_file $tissue_name $sig_thresh $variant_annotation_enrichment_file


    pip_thresh="0.5"
    sig_thresh="0.15"
    processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}"PIP_"${pip_thresh}"_fine_mapped_eqtl_results.txt"
    borzoi_results_stem=$organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty_"
    variant_annotation_enrichment_file=$variant_annotation_enrichment_dir"variant_enrichment_in_fine_mapped_eqtls_"${tissue_name}"_"${sig_thresh}"_"${pip_thresh}"_summary.txt"
    sbatch borzoi_variant_enrichment_in_fm_eqtls.sh $borzoi_results_stem $processed_fm_eqtl_output_file $eqtl_sumstats_dir $gtex_eqtl_tissue_names_file $tissue_name $sig_thresh $variant_annotation_enrichment_file
done
fi


if false; then
sig_thresh="0.05"
tissue_overlap_output_file=$variant_annotation_enrichment_dir"tisue_overlaps_"${sig_thresh}".txt"
sbatch compute_tissue_overap_statistics.sh $sig_thresh $tissue_overlap_output_file $gtex_tissue_names_file $organized_borzoi_gtex_predictions

sig_thresh="0.1"
tissue_overlap_output_file=$variant_annotation_enrichment_dir"tisue_overlaps_"${sig_thresh}".txt"
sbatch compute_tissue_overap_statistics.sh $sig_thresh $tissue_overlap_output_file $gtex_tissue_names_file $organized_borzoi_gtex_predictions
fi


if false; then
source ~/.bashrc
conda activate borzoi
python visualize_variant_anno_enrichment.py $visualize_variant_anno_enrich_dir $variant_annotation_enrichment_dir $gtex_tissue_names_file $non_redundant_gwas_traits_file $gwas_sldsc_results_dir
fi

if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_variant_enrichments.R $visualize_variant_anno_enrich_dir $variant_annotation_enrichment_dir $gtex_eqtl_tissue_names_file
fi




###################
# Trait h2 enrichments by borzoi results
tissue_name="Muscle_Skeletal"
borzoi_res_file=$organized_borzoi_gtex_predictions${tissue_name}"_borzoi_estimates_w_uncertainty.txt"
if false; then
sbatch create_sldsc_annotation_files_for_single_tissue.sh $tissue_name $sldsc_processed_anno_dir $borzoi_res_file $ldsc_snp_annotation_dir $genotype_1000G_plink_stem ${ldsc_code_dir} $ldsc_weights_dir
fi

if false; then
source ~/.bashrc
conda activate ldsc
trait_name="UKB_460K.disease_AID_ALL"
gwas_file_name=${gwas_sumstats_dir}${trait_name}".sumstats"
python ${ldsc_code_dir}ldsc.py --h2 ${gwas_file_name} --ref-ld-chr ${ldsc_snp_annotation_dir}"baselineLD." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${genotype_1000G_plink_stem} --out ${sldsc_h2_results_dir}${trait_name}"_sldsc_testing"
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


