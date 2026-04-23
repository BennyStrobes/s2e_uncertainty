#################
# Input data
#################


# Directory containing results of borzoi runs
borzoi_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/"

# Directory containing borzoi gtex target indices and names
borzoi_gtex_target_names_file=${borzoi_results_dir}"targets_gtex_only_ordered.txt"
borzoi_gtex_independent_target_names_file=${borzoi_results_dir}"targets_gtex_only_independent_ordered.txt"

# Directory containing genotype data
processed_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/plink_processed_genotype/"

# Directory of expression data
gtex_expr_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/residualized_expression/"

# Directory containing eQTL summary statistics
eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/eqtl_results/"

# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"



#################
# Output directories
#################
# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/eqtl_borzoi_correlations/"

borzoi_output_dir=${output_root}"processed_borzoi/"
ld_corr_results_output_dir=${output_root}"ld_corr_results/"

borzoi_based_prior_output_dir=${output_root}"borzoi_based_prior/"

bayesian_ld_corr_processed_data_dir=${output_root}"bayesian_ld_corr_processed_data/"
bayesian_ld_corr_results_dir=${output_root}"bayesian_ld_corr_results/"



visualization_dir=${output_root}"visualization/"



#################
# Code
#################

#################
# Preprocess data to get into correct format for LD corr
#################

#####
# 1. eQTLs have already been processed

#####
# 2. Borzoi effects
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	sbatch preprocess_borzoi_data_for_ld_corr.sh $gtex_v10_pc_genes_gtf $borzoi_results_dir $borzoi_target_index $gtex_tissue $target_identifier $borzoi_output_dir
done
fi


#####
# 3. Annotation effects
if false; then
anno_methods=("borzoi_magnitude_bins" "dist_to_tss_bins" "strand_dist_to_tss_bins" "af_bins" "borzoi_magnitude_binsXaf_bins" "borzoi_magnitude_binsXdist_to_tss_bins")
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	for anno_method in "${anno_methods[@]}"; do
		borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
		sbatch annotate_variant_gene_pairs.sh $borzoi_effect_file $anno_method $borzoi_annotation_file $eqtl_sumstats_file $target_tissue
	done
done
fi








#################
# Run LD-corr
#################
if false; then
anno_methods=("borzoi_magnitude_bins" "dist_to_tss_bins" "strand_dist_to_tss_bins" "af_bins" "borzoi_magnitude_binsXaf_bins" "borzoi_magnitude_binsXdist_to_tss_bins")
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	for anno_method in "${anno_methods[@]}"; do
		borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
		ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_"${target_tissue}"_"${target_sample}"_"${anno_method}
		sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_sumstats_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $ld_corr_output_stem
	done
done
fi








#################
# Run borzoi prior inference
#################
if false; then
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="gaussian"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}


	sh run_borzoi_based_prior_inference.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem
fi

















####################
# OLD
#####################

target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"

if false; then
eqtl_effects_file=$eqtl_ss_output_dir${target_tissue}"_eqtl_sumstats.txt.gz"
borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
genotype_stem=$genotype_output_dir"gtex_v9_eqtl_chr"

ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem

target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="dist_to_tss_bins"

eqtl_effects_file=$eqtl_ss_output_dir${target_tissue}"_eqtl_sumstats.txt.gz"
borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
genotype_stem=$genotype_output_dir"gtex_v9_eqtl_chr"

ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem


target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="af_bins"

eqtl_effects_file=$eqtl_ss_output_dir${target_tissue}"_eqtl_sumstats.txt.gz"
borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
genotype_stem=$genotype_output_dir"gtex_v9_eqtl_chr"

ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem

fi



target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"

eqtl_effects_file=$eqtl_ss_output_dir${target_tissue}"_eqtl_sumstats.txt.gz"
borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
genotype_stem=$genotype_output_dir"gtex_v9_eqtl_chr"

if false; then
chi_sq_thresh="80"
ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${chi_sq_thresh}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem $chi_sq_thresh

chi_sq_thresh="200"
ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${chi_sq_thresh}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem $chi_sq_thresh

chi_sq_thresh="400"
ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${chi_sq_thresh}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem $chi_sq_thresh

chi_sq_thresh="1000"
ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${chi_sq_thresh}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem $chi_sq_thresh
fi


if false; then
anno_method="af_bins"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
chi_sq_thresh="80"
ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${chi_sq_thresh}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem $chi_sq_thresh

anno_method="dist_to_tss_bins"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
chi_sq_thresh="80"
ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_std_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${chi_sq_thresh}
sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_effects_file $borzoi_annotation_file $genotype_stem $ld_corr_output_stem $chi_sq_thresh
fi


#################
# Run LD-corr compete
#################
#####
# 2. Create Borzoi files for this analysis
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	sbatch preprocess_borzoi_data_for_ld_corr_compete.sh $gtex_v10_pc_genes_gtf $borzoi_results_dir $borzoi_target_index $gtex_tissue $target_identifier $borzoi_output_dir

done
fi



# Create file summarizing all Borzoi files (one file/tissue for each line)
borzoi_tracks_file=${borzoi_output_dir}"borzoi_tracks.txt"
if false; then
echo -e "track_name\tborzoi_effect_file" > $borzoi_tracks_file
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	echo -e "${gtex_tissue}\t${borzoi_output_dir}${gtex_tissue}_${target_identifier}_borzoi_effects_compete_version.txt.gz" >> $borzoi_tracks_file
done
fi

# Run LD-Corr compete analysis
if false; then
model="joint" # joint vs independent
thresh="0.0"
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	eqtl_effects_file=$eqtl_ss_output_dir${gtex_tissue}"_eqtl_sumstats.txt.gz"
	ld_corr_compete_output_stem=${ld_corr_results_output_dir}"ld_corr_compete_results_"$model"_"${thresh}"_"${gtex_tissue}"_chisq_thresh"
	sbatch run_ld_corr_compete.sh $borzoi_tracks_file $eqtl_effects_file $genotype_stem $ld_corr_compete_output_stem $model $thresh
done
fi
if false; then

model="independent" # joint vs independent
thresh="0.05"
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	eqtl_effects_file=$eqtl_ss_output_dir${gtex_tissue}"_eqtl_sumstats.txt.gz"
	ld_corr_compete_output_stem=${ld_corr_results_output_dir}"ld_corr_compete_results_"$model"_"${thresh}"_"${gtex_tissue}
	sbatch run_ld_corr_compete.sh $borzoi_tracks_file $eqtl_effects_file $genotype_stem $ld_corr_compete_output_stem $model $thresh
done
fi



#################
# Bayesian LD corr
#################
bayes_input_data_stem=$bayesian_ld_corr_processed_data_dir"ten_independent_tissues_data_"
if false; then
sh preprocess_data_for_bayesian_ld_corr.sh $borzoi_gtex_independent_target_names_file $gtex_v10_pc_genes_gtf $eqtl_sumstats_dir $borzoi_results_dir $genotype_stem $bayes_input_data_stem
fi

# LD CORR compete
if false; then
bayes_input_data_summary_file=$bayesian_ld_corr_processed_data_dir"ten_independent_tissues_data__per_gene_summary.txt"
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description tissue_name; do

	prior_version="non_informative_priors"
	bayes_ld_corr_output_stem=$bayesian_ld_corr_results_dir"bayesian_ld_corr_compete_"${tissue_name}"_"${prior_version}
	sbatch run_bayes_ld_corr_compete.sh $bayes_input_data_summary_file $tissue_name $borzoi_gtex_independent_target_names_file $bayes_ld_corr_output_stem $prior_version

	prior_version="ard_prior_means"
	bayes_ld_corr_output_stem=$bayesian_ld_corr_results_dir"bayesian_ld_corr_compete_"${tissue_name}"_"${prior_version}
	sbatch run_bayes_ld_corr_compete.sh $bayes_input_data_summary_file $tissue_name $borzoi_gtex_independent_target_names_file $bayes_ld_corr_output_stem $prior_version
done
fi




# LD corr stratified
if false; then
bayes_input_data_summary_file=$bayesian_ld_corr_processed_data_dir"ten_independent_tissues_data__per_gene_summary.txt"
tissue_name="Adipose_Subcutaneous"
	bayes_ld_corr_output_stem=$bayesian_ld_corr_results_dir"bayesian_ld_corr_magnitude_stratified_"${tissue_name}
	sh run_bayes_ld_corr_magnitude_stratified.sh $bayes_input_data_summary_file $tissue_name $borzoi_gtex_independent_target_names_file $bayes_ld_corr_output_stem
fi






#################
# Visualize results
#################
if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_correlation_results.R ${ld_corr_results_output_dir} $visualization_dir
fi

















