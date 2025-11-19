##################
# Input data
##################
# Downloaded from (on 10/20/25):
# curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz  
# gunzip hg19.fa.gz
fasta_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/hg19.fa"
# Downloaded from (on 10/20/25):
# https://github.com/ni-lab/finetuning-enformer/tree/main/process_geuvadis_data/log_tpm/corrected_log_tpm.annot.csv.gz
expr_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/GEUVADIS/corrected_log_tpm.annot.csv.gz"
# Downloaded from https://github.com/mostafavilabuw/SAGEnet/tree/main/input_data on 10/20/25
tss_data_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/gene-ids-and-positions.tsv"
# Downloaded from https://github.com/mostafavilabuw/SAGEnet/tree/main/input_data on 10/20/25
protein_coding_gene_list="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/protein_coding_genes.csv"

# Gtex fine-mapped eqtl results
gtex_fine_mapped_eqtl_file="/lab-share/CHIP-Strober-e2/Public/GTEx/fine_mapping/GTEx_49tissues_release1.tsv"

# Downloaded from (on 10/20/25):
# curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz  
# gunzip hg38.fa.gz
hg38_fasta_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/hg38.fa"

gtex_wb_tpm_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/GTEX/gene_tpm_2017-06-05_v8_whole_blood.gct.gz"
gtex_cortex_tpm_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/GTEX/gene_tpm_2017-06-05_v8_brain_cortex.gct.gz"


gtex_cortex_log_tpm_from_anna="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/expression_from_anna/GTEx_Brain-Cortex_covariate_adjusted_log_tpm_from_anna.tsv"
rosmap_cortex_log_tpm_from_anna="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/expression_from_anna/Rosmap_vcf_match_covariate_adjusted_log_tpm.tsv"


##################
# Output data
##################
# Directory to keep track of model fits
model_training_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/model_training/"
# directory containing evaluation results
model_evaluation_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/model_evaluation/"

processed_expression_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/processed_expression/"

processed_fm_eqtl_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/processed_fm_eqtl_data_dir/"

delta_score_fm_eqtl_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/delta_score_fm_eqtl/"


##################
# Load in environment to run SAGE
##################
if false; then
conda activate SAGEnet
fi



##################
# Run analysis
##################
# Train model based on rsage notebook
if false; then
sh train_rsage_based_on_notebook.sh $fasta_file $expr_file $tss_data_file $protein_coding_gene_list $model_training_dir
fi

# Evaluate model based on rsage notebook
if false; then
sh evaluate_rsage_based_on_notebook.sh $fasta_file $expr_file $tss_data_file $protein_coding_gene_list $model_training_dir $model_evaluation_dir
fi

# Process gtex expression data
if false; then
tpm_thresh="0.1"
max_prop_sample_missing="0.2"
gtex_wb_processed_expression_file_stem=${processed_expression_dir}"GTEx_whole_blood_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}
sh process_gtex_expression_data.sh $gtex_wb_tpm_file $gtex_wb_processed_expression_file_stem $tpm_thresh $max_prop_sample_missing

tpm_thresh="0.1"
max_prop_sample_missing="0.5"
gtex_wb_processed_expression_file_stem=${processed_expression_dir}"GTEx_whole_blood_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}
sh process_gtex_expression_data.sh $gtex_wb_tpm_file $gtex_wb_processed_expression_file_stem $tpm_thresh $max_prop_sample_missing

tpm_thresh="0.0"
max_prop_sample_missing="0.5"
gtex_wb_processed_expression_file_stem=${processed_expression_dir}"GTEx_whole_blood_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}
sh process_gtex_expression_data.sh $gtex_wb_tpm_file $gtex_wb_processed_expression_file_stem $tpm_thresh $max_prop_sample_missing
fi



# GTEx
if false; then
tpm_thresh="0.0"
max_prop_sample_missing="0.5"
gtex_wb_processed_expression_file=${processed_expression_dir}"GTEx_whole_blood_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
output_stem="GTEx_Whole_Blood_"${tpm_thresh}"_"${max_prop_sample_missing}
sbatch train_rsage_on_gtex.sh $hg38_fasta_file $gtex_wb_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir $output_stem

tpm_thresh="0.1"
max_prop_sample_missing="0.2"
gtex_wb_processed_expression_file=${processed_expression_dir}"GTEx_whole_blood_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
output_stem="GTEx_Whole_Blood_"${tpm_thresh}"_"${max_prop_sample_missing}
sbatch train_rsage_on_gtex.sh $hg38_fasta_file $gtex_wb_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir $output_stem
fi



# Evaluate model based on rsage notebook
tpm_thresh="0.0"
max_prop_sample_missing="0.5"
gtex_wb_processed_expression_file=${processed_expression_dir}"GTEx_whole_blood_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
output_stem="GTEx_Whole_Blood_"${tpm_thresh}"_"${max_prop_sample_missing}
if false; then
sh evaluate_rsage_on_gtex.sh $hg38_fasta_file $gtex_wb_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir $model_evaluation_dir $output_stem
fi

# Parse eQTL data
tissue_name="Whole_Blood"
pip_threshold="0.9"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}$tissue_name"_PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"
if false; then
sh parse_eqtl_data.sh $gtex_fine_mapped_eqtl_file $tissue_name $pip_threshold $processed_fm_eqtl_output_file
fi

# Get delta scores for fine-mapped eqtl variants
if false; then
sh get_delta_scores_for_fm_eqtl_variants.sh $processed_fm_eqtl_output_file $hg38_fasta_file $gtex_wb_tpm_file $tss_data_file $protein_coding_gene_list $model_training_dir'epoch=15-step=230608.ckpt' $delta_score_fm_eqtl_dir"Whole_Blood_0.9_summary.txt"
sh get_delta_scores_for_fm_eqtl_variants.sh $processed_fm_eqtl_output_file $hg38_fasta_file $gtex_wb_tpm_file $tss_data_file $protein_coding_gene_list $model_training_dir'epoch=13-step=156366.ckpt' $delta_score_fm_eqtl_dir"Whole_Blood_0.9_summary.txt"
fi


###########
# CORTEX
###########
# Process gtex expression data
if false; then
tpm_thresh="0.1"
max_prop_sample_missing="0.2"
gtex_cortex_processed_expression_file_stem=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}
sh process_gtex_expression_data.sh $gtex_cortex_tpm_file $gtex_cortex_processed_expression_file_stem $tpm_thresh $max_prop_sample_missing

tpm_thresh="0.1"
max_prop_sample_missing="0.5"
gtex_cortex_processed_expression_file_stem=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}
sh process_gtex_expression_data.sh $gtex_cortex_tpm_file $gtex_cortex_processed_expression_file_stem $tpm_thresh $max_prop_sample_missing

tpm_thresh="0.0"
max_prop_sample_missing="0.5"
gtex_cortex_processed_expression_file_stem=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}
sh process_gtex_expression_data.sh $gtex_cortex_tpm_file $gtex_cortex_processed_expression_file_stem $tpm_thresh $max_prop_sample_missing
fi


if false; then
tpm_thresh="0.0"
max_prop_sample_missing="0.5"
gtex_cortex_processed_expression_file=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
output_stem="GTEx_cortex_"${tpm_thresh}"_"${max_prop_sample_missing}
sbatch train_rsage_on_gtex.sh $hg38_fasta_file $gtex_cortex_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir $output_stem



tpm_thresh="0.1"
max_prop_sample_missing="0.2"
gtex_cortex_processed_expression_file=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
output_stem="GTEx_cortex_"${tpm_thresh}"_"${max_prop_sample_missing}
sbatch train_rsage_on_gtex.sh $hg38_fasta_file $gtex_cortex_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir $output_stem
fi


# Evaluate model based on rsage notebook
tpm_thresh="0.0"
max_prop_sample_missing="0.5"
gtex_cortex_processed_expression_file=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
output_stem="GTEx_cortex_"${tpm_thresh}"_"${max_prop_sample_missing}
if false; then
sh evaluate_rsage_on_gtex.sh $hg38_fasta_file $gtex_cortex_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir $model_evaluation_dir $output_stem
fi

tpm_thresh="0.1"
max_prop_sample_missing="0.2"
gtex_cortex_processed_expression_file=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
output_stem="GTEx_cortex_"${tpm_thresh}"_"${max_prop_sample_missing}
if false; then
sh evaluate_rsage_on_gtex.sh $hg38_fasta_file $gtex_cortex_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir $model_evaluation_dir $output_stem
fi




# Parse eQTL data
tissue_name="Brain_Cortex"
pip_threshold="0.9"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}$tissue_name"_PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"
if false; then
sh parse_eqtl_data.sh $gtex_fine_mapped_eqtl_file $tissue_name $pip_threshold $processed_fm_eqtl_output_file
fi


# Get delta scores for fine-mapped eqtl variants

tpm_thresh="0.0"
max_prop_sample_missing="0.5"
gtex_processed_expression_file=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
bootstrap_iter="1"
output_stem="GTEx_cortex_"${tpm_thresh}"_"${max_prop_sample_missing}"_bs_"${bootstrap_iter}
if false; then
sh get_delta_scores_for_fm_eqtl_variants.sh $processed_fm_eqtl_output_file $hg38_fasta_file $gtex_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir${output_stem}"/" $delta_score_fm_eqtl_dir"Cortex_0.9_summary.txt"
fi



##################
# Try bootstrapping
##################

if false; then
tpm_thresh="0.0"
max_prop_sample_missing="0.5"
gtex_cortex_processed_expression_file=${processed_expression_dir}"GTEx_cortex_log_tpm_"${tpm_thresh}"_"${max_prop_sample_missing}".txt"
for bootstrap_iter in {1..20}; do
output_stem="GTEx_cortex_"${tpm_thresh}"_"${max_prop_sample_missing}"_bs_"${bootstrap_iter}
sbatch train_rsage_on_gtex_bootstrapping.sh $hg38_fasta_file $gtex_cortex_processed_expression_file $tss_data_file $protein_coding_gene_list $model_training_dir $output_stem $bootstrap_iter
done
fi



# Parse eQTL data
tissue_name="Brain_Cortex"
pip_threshold="0.9"
processed_fm_eqtl_output_file=${processed_fm_eqtl_data_dir}$tissue_name"_PIP_"${pip_threshold}"_fine_mapped_eqtl_results.txt"

tpm_thresh="0.0"
max_prop_sample_missing="0.5"
n_bootstraps="20"
model_output_stem="GTEx_cortex_"${tpm_thresh}"_"${max_prop_sample_missing}"_bs_"
if false; then
sh get_bs_delta_scores_for_fm_eqtl_variants.sh $processed_fm_eqtl_output_file $hg38_fasta_file $tss_data_file $protein_coding_gene_list $model_training_dir${model_output_stem} $delta_score_fm_eqtl_dir"Cortex_0.9_summary.txt" $n_bootstraps
fi




##############################
# Use anna's downloaded files
###############################
if false; then
# USE GTEx expression cortex from anna for training
output_stem="GTEx_cortex_from_anna"
sbatch train_rsage_on_gtex.sh $hg38_fasta_file $gtex_cortex_log_tpm_from_anna $tss_data_file $protein_coding_gene_list $model_training_dir $output_stem


# USE Rosmap expression cortex from anna for training
output_stem="Rosmap_cortex_from_anna"
sbatch train_rsage_on_gtex.sh $hg38_fasta_file $rosmap_cortex_log_tpm_from_anna $tss_data_file $protein_coding_gene_list $model_training_dir $output_stem

fi



# Evaluate model based on rsage notebook
if false; then
output_stem="GTEx_cortex_from_anna"
sh evaluate_rsage_on_gtex.sh $hg38_fasta_file $gtex_cortex_log_tpm_from_anna $tss_data_file $protein_coding_gene_list $model_training_dir $model_evaluation_dir $output_stem

output_stem="Rosmap_cortex_from_anna"
sh evaluate_rsage_on_gtex.sh $hg38_fasta_file $rosmap_cortex_log_tpm_from_anna $tss_data_file $protein_coding_gene_list $model_training_dir $model_evaluation_dir $output_stem
fi


if false; then
output_stem="GTEx_cortex_from_anna"
sh get_delta_scores_for_fm_eqtl_variants.sh $processed_fm_eqtl_output_file $hg38_fasta_file $gtex_cortex_log_tpm_from_anna $tss_data_file $protein_coding_gene_list $model_training_dir${output_stem}"/" $delta_score_fm_eqtl_dir${output_stem}".txt"


output_stem="Rosmap_cortex_from_anna"
sh get_delta_scores_for_fm_eqtl_variants.sh $processed_fm_eqtl_output_file $hg38_fasta_file $rosmap_cortex_log_tpm_from_anna $tss_data_file $protein_coding_gene_list $model_training_dir${output_stem}"/" $delta_score_fm_eqtl_dir${output_stem}".txt"
fi
























