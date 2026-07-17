#################
# Input data
#################


# Directory containing results of borzoi runs
borzoi_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/"

# Directory containing borzoi gtex target indices and names
borzoi_gtex_target_names_file=${borzoi_results_dir}"targets_gtex_only_ordered.txt"
borzoi_gtex_unique_target_names_file=${borzoi_results_dir}"targets_gtex_eqtl_only_unique_ordered.txt"

borzoi_gtex_independent_target_names_file=${borzoi_results_dir}"targets_gtex_only_independent_ordered.txt"

borzoi_non_gtex_target_names_file=${borzoi_results_dir}"targets_interesting_non_gtex_ordered.txt"



# Directory containing genotype data
processed_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/plink_processed_genotype/"

# Directory of expression data
gtex_expr_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/residualized_expression/"

# Directory containing eQTL summary statistics
eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/eqtl_results/"

# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"

# GWAS sumstats dir
gwas_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/UKBB_all_snps_sumstats/data/"

simulation_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/simulations/corr_sim_inference_results/"

baselineLD_anno_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

gtex_fm_file="/lab-share/CHIP-Strober-e2/Public/GTEx/fine_mapping/v8/GTEx_49tissues_release1.tsv"

#################
# Output directories
#################
# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/eqtl_borzoi_correlations/"

borzoi_output_dir=${output_root}"processed_borzoi/"

bootstrapped_cross_tissue_gene_sets_dir=${output_root}"bootstrapped_gene_sets/"

ld_corr_results_output_dir=${output_root}"ld_corr_results/"

cross_tissue_ld_corr_results_output_dir=${output_root}"cross_tissue_ld_corr_results/"

borzoi_based_prior_output_dir=${output_root}"borzoi_based_prior/"

expr_pred_output_dir=${output_root}"expression_prediction/"

visualize_expression_pred=${output_root}"visualize_expression_pred/"

bayesian_ld_corr_processed_data_dir=${output_root}"bayesian_ld_corr_processed_data/"
bayesian_ld_corr_results_dir=${output_root}"bayesian_ld_corr_results/"


hybrid_expression_prediction_dir=${output_root}"hybrid_expression_prediction/"

new_hybrid_expression_prediction_dir=${output_root}"new_hybrid_expression_prediction/"

hybrid_coloc_dir=${output_root}"coloc_hybrid/"


twas_dir=${output_root}"twas/"

visualize_ld_corr_results_dir=${output_root}"visualize_ld_corr/"

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
tail -n +2 "$borzoi_gtex_unique_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	sbatch preprocess_borzoi_data_for_ld_corr.sh $gtex_v10_pc_genes_gtf $borzoi_results_dir $borzoi_target_index $gtex_tissue $target_identifier $borzoi_output_dir
done
fi

if false; then
tail -n +2 "$borzoi_non_gtex_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	sbatch preprocess_borzoi_data_for_ld_corr.sh $gtex_v10_pc_genes_gtf $borzoi_results_dir $borzoi_target_index $gtex_tissue $target_identifier $borzoi_output_dir
done
fi



#####
# 3. Annotation effects
anno_methods=("borzoi_magnitude_bins" "dist_to_tss_bins" "strand_dist_to_tss_bins" "af_bins" "borzoi_magnitude_binsXaf_bins" "borzoi_magnitude_binsXdist_to_tss_bins" "borzoi_effect_size_bins" "borzoi_finer_effect_size_bins")
anno_methods=("intercept")

anno_methods=("baselineLD_Conserved_Mammal_phastCons46way" "baselineLD_TSS_Hoffman" "baselineLD_DHS_Trynka" "baselineLD_H3K27ac_Hnisz" "baselineLD_Enhancer_Hoffman")
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	for anno_method in "${anno_methods[@]}"; do
		borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
		sbatch annotate_variant_gene_pairs.sh $borzoi_effect_file $anno_method $borzoi_annotation_file $eqtl_sumstats_file $target_tissue $baselineLD_anno_dir
	done
done
fi


#################
# Generate cross tissue gene sets (and bootstrapped gene sets)
#################
if false; then
sbatch generate_cross_tissue_bootstrapped_gene_sets.sh ${eqtl_sumstats_dir} ${borzoi_output_dir} ${borzoi_gtex_unique_target_names_file} ${bootstrapped_cross_tissue_gene_sets_dir}
fi


#################
# Run LD-corr
#################


target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"

eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"
if false; then
anno_method="borzoi_magnitude_bins"
		borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
		ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_newer"
		sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_sumstats_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file ${bootstrapped_cross_tissue_gene_sets_dir}"cross_tissue_gene_set_bootstrap_" $ld_corr_output_stem
fi





anno_methods=("borzoi_magnitude_bins" "dist_to_tss_bins" "strand_dist_to_tss_bins" "af_bins" "borzoi_magnitude_binsXaf_bins" "borzoi_magnitude_binsXdist_to_tss_bins")
anno_methods=("intercept")

if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	echo $target_tissue" "$target_sample
	for anno_method in "${anno_methods[@]}"; do
		borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
		ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_"${target_tissue}"_"${target_sample}"_"${anno_method}
		sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_sumstats_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $ld_corr_output_stem
	done
done
fi
if false; then
anno_methods=("baselineLD_Conserved_Mammal_phastCons46way" "baselineLD_TSS_Hoffman" "baselineLD_DHS_Trynka" "baselineLD_H3K27ac_Hnisz" "baselineLD_Enhancer_Hoffman")
anno_methods=("borzoi_magnitude_bins" "dist_to_tss_bins")
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	echo $target_tissue" "$target_sample
	for anno_method in "${anno_methods[@]}"; do
		borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
		ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_"${target_tissue}"_"${target_sample}"_"${anno_method}
		sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_sumstats_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $ld_corr_output_stem
	done
fi


target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"

eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

anno_method="borzoi_magnitude_bins"

anno_methods=("intercept" "borzoi_magnitude_bins" "dist_to_tss_bins" "strand_dist_to_tss_bins" "af_bins" "borzoi_magnitude_binsXaf_bins" "borzoi_magnitude_binsXdist_to_tss_bins" "baselineLD_Conserved_Mammal_phastCons46way" "baselineLD_TSS_Hoffman" "baselineLD_DHS_Trynka" "baselineLD_H3K27ac_Hnisz" "baselineLD_Enhancer_Hoffman")
if false; then
for anno_method in "${anno_methods[@]}"; do
		borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
		ld_corr_output_stem=${ld_corr_results_output_dir}"ld_corr_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_newer"
		sbatch run_ld_corr.sh $borzoi_effect_file $eqtl_sumstats_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $ld_corr_output_stem
done
fi








#################
# Run Cross trait LD-corr
#################
if false; then
anno_methods="borzoi_magnitude_bins intercept"
genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"

t1_index=0
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS="$(printf '\t')" read -r t1_orig_target_index t1_borzoi_target_index t1_target_sample t1_target_description t1_target_tissue; do
	t1_eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${t1_target_tissue}"_sumstats.txt.gz"
	t1_borzoi_effect_file=${borzoi_output_dir}${t1_target_tissue}"_"${t1_target_sample}"_borzoi_effects.txt.gz"
	t1_genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${t1_target_tissue}"_expression_samples.txt"

	t2_index=0
	tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS="$(printf '\t')" read -r t2_orig_target_index t2_borzoi_target_index t2_target_sample t2_target_description t2_target_tissue; do
		if [ "$t2_index" -le "$t1_index" ]; then
			t2_index=$((t2_index + 1))
			continue
		fi

		t2_eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${t2_target_tissue}"_sumstats.txt.gz"
		t2_borzoi_effect_file=${borzoi_output_dir}${t2_target_tissue}"_"${t2_target_sample}"_borzoi_effects.txt.gz"
		t2_genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${t2_target_tissue}"_expression_samples.txt"

		echo $t1_target_tissue" "$t1_target_sample" "$t2_target_tissue" "$t2_target_sample
		for anno_method in $anno_methods; do
			t1_borzoi_annotation_file=${borzoi_output_dir}${t1_target_tissue}"_"${t1_target_sample}"_"${anno_method}"_annotations.txt.gz"
			t2_borzoi_annotation_file=${borzoi_output_dir}${t2_target_tissue}"_"${t2_target_sample}"_"${anno_method}"_annotations.txt.gz"
			cross_tissue_ld_corr_output_stem=${cross_tissue_ld_corr_results_output_dir}"ld_corr_results_"${t1_target_tissue}"_"${t1_target_sample}"_"${t2_target_tissue}"_"${t2_target_sample}"_"${anno_method}

			sbatch run_cross_tissue_ld_corr.sh $t1_borzoi_effect_file $t2_borzoi_effect_file $t1_eqtl_sumstats_file $t2_eqtl_sumstats_file $t1_borzoi_annotation_file $t2_borzoi_annotation_file $genotype_stem $t1_genotype_sample_mapping_file $t2_genotype_sample_mapping_file $cross_tissue_ld_corr_output_stem
		done
		t2_index=$((t2_index + 1))
	done
	t1_index=$((t1_index + 1))
done
fi


#################
# Visualize results
#################
if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_correlation_results.R ${ld_corr_results_output_dir} $simulation_results_dir $cross_tissue_ld_corr_results_output_dir $visualize_ld_corr_results_dir
fi





#################
# Run borzoi prior inference
#################
if false; then
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="ashr"

standardize_geno="False"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}

	sbatch run_borzoi_based_prior_inference.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno

standardize_geno="True"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}

	sbatch run_borzoi_based_prior_inference.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
fi

if false; then
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="gaussian"

	standardize_geno="False"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}


	sbatch run_borzoi_based_prior_inference.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno


target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="gaussian"

	standardize_geno="True"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}

	sbatch run_borzoi_based_prior_inference.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
fi

if false; then
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="borzoi_scaled_variance"

	standardize_geno="True"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}

	sh run_borzoi_based_prior_inference.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
fi


if false; then
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_effect_size_bins"
distribution="uniform_prior_grid"

	standardize_geno="False"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"grid2alpha.1_"${standardize_geno}

	sbatch run_borzoi_based_prior_inference.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $standardize_geno
fi
if false; then
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_effect_sizes"


	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"

	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"




	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_two_step_gradient_based_"${target_tissue}"_"${target_sample}"_"${anno_method}

	sbatch run_borzoi_based_prior_inference_gradient_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem
fi


target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_effect_sizes"
model_effect_scale="allelic"

if false; then
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"

	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"





	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_squared_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem

	


target_tissue="Whole_Blood"
target_sample="GTEX-1LB8K-0005-SM-DIPED.1"
anno_method="borzoi_effect_sizes"
model_effect_scale="allelic"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"

	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_squared_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem


target_tissue="Liver"
target_sample="GTEX-11EQ9-0526-SM-5A5JZ.1"
anno_method="borzoi_effect_sizes"
model_effect_scale="allelic"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"

	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_squared_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem

target_tissue="Artery_Aorta"
target_sample="GTEX-1JK1U-0426-SM-CYPSP.1"
anno_method="borzoi_effect_sizes"
model_effect_scale="allelic"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"

	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_squared_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem

target_tissue="Adipose_Subcutaneous"
target_sample="GTEX-132QS-2526-SM-62LFJ.1"
anno_method="borzoi_effect_sizes"
model_effect_scale="allelic"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"

	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_squared_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem

target_tissue="Lung"
target_sample="GTEX-1399S-1726-SM-5L3DI.1"
anno_method="borzoi_effect_sizes"
model_effect_scale="allelic"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"

	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
	sbatch run_borzoi_based_prior_inference_ldscore_grid_squared_based.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $eqtl_sumstats_file ${anno_method} $borzoi_based_prior_output_stem

fi





if false; then
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="gaussian"
standardize_geno="False"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}

	expr_prediction_output_file=${expr_pred_output_dir}"expr_pred_summary_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_full.txt"
	sbatch predict_cross_individual_expression.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno


target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="gaussian"
standardize_geno="True"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}

	expr_prediction_output_file=${expr_pred_output_dir}"expr_pred_summary_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_full.txt"
	sbatch predict_cross_individual_expression.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno
fi




########### Predict expression with Whole Blood LDSC squared-grid prior
if false; then
target_tissue="Whole_Blood"
target_sample="GTEX-1LB8K-0005-SM-DIPED.1"
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}

	expr_prediction_output_file=${expr_pred_output_dir}"expr_pred_summary_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt"
	sh predict_cross_individual_expression.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno
fi

if false; then
ref_tissue="Whole_Blood"
ref_sample="GTEX-1LB8K-0005-SM-DIPED.1"
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}

	expr_prediction_output_file=${expr_pred_output_dir}"expr_pred_summary_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt"
	sbatch predict_cross_individual_expression.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno
done
fi





if false; then
ref_tissue="Muscle_Skeletal"
ref_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"

	target_tissue="Whole_Blood"
	target_sample="GTEX-1LB8K-0005-SM-DIPED.1"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}

	expr_prediction_output_file=${expr_pred_output_dir}"expr_pred_summary_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt"
	sbatch predict_cross_individual_expression.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno
fi
if false; then
ref_tissue="Muscle_Skeletal"
ref_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid"
standardize_geno="False"

	target_tissue="Whole_Blood"
	target_sample="GTEX-1LB8K-0005-SM-DIPED.1"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}

	expr_prediction_output_file=${expr_pred_output_dir}"expr_pred_summary_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt"
	sbatch predict_cross_individual_expression.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno
fi


#################
# Fine-mapped variant prediction
if false; then
target_tissue="Whole_Blood"
target_sample="GTEX-1LB8K-0005-SM-DIPED.1"
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"


	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}

	fine_mapping_method="SUSIE"
	pip_threshold="0.9"
	finemapped_borzoi_prior_output_file=${borzoi_based_prior_output_dir}"fine_mapped_eqtl_borzoi_prior_comparison_"${target_tissue}"_"${target_sample}"_"${fine_mapping_method}"_pip_"${pip_threshold}".txt"
	sh compare_finemapped_eqtls_to_borzoi_prior.sh $gtex_fm_file $target_tissue $fine_mapping_method $pip_threshold $borzoi_effect_file $borzoi_based_prior_output_stem $finemapped_borzoi_prior_output_file
fi









#################
# Marginal effect prediction
if false; then
ref_tissue="Whole_Blood"
ref_sample="GTEX-1LB8K-0005-SM-DIPED.1"
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}

	marginal_effect_prediction_output_file=${expr_pred_output_dir}"marginal_eqtl_effect_prediction_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt.gz"
	ld_prune_r_sq_threshold="0.2"
	ld_prune_random_seed="1"
	ld_pruned_marginal_effect_prediction_output_file=${expr_pred_output_dir}"marginal_eqtl_effect_prediction_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_ld_pruned_r2_"${ld_prune_r_sq_threshold}"_seed_"${ld_prune_random_seed}".txt.gz"
	sbatch predict_marginal_eqtl_effects_from_borzoi_prior.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $borzoi_based_prior_output_stem $marginal_effect_prediction_output_file $ld_pruned_marginal_effect_prediction_output_file $ld_prune_r_sq_threshold $ld_prune_random_seed
done
fi

ref_tissue="Muscle_Skeletal"
ref_sample="GTEX-13QJ3-0726-SM-5SI68.1"

target_tissue="Whole_Blood"
target_sample="GTEX-1LB8K-0005-SM-DIPED.1"

anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"
if false; then

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}

	marginal_effect_prediction_output_file=${expr_pred_output_dir}"marginal_eqtl_effect_prediction_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt.gz"
	ld_prune_r_sq_threshold="0.2"
	ld_prune_random_seed="1"
	ld_pruned_marginal_effect_prediction_output_file=${expr_pred_output_dir}"marginal_eqtl_effect_prediction_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_ld_pruned_r2_"${ld_prune_r_sq_threshold}"_seed_"${ld_prune_random_seed}".txt.gz"
	sbatch predict_marginal_eqtl_effects_from_borzoi_prior.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $borzoi_based_prior_output_stem $marginal_effect_prediction_output_file $ld_pruned_marginal_effect_prediction_output_file $ld_prune_r_sq_threshold $ld_prune_random_seed


anno_method="borzoi_effect_sizes"
distribution="ldscore_grid"
standardize_geno="False"

	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}
	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}



	marginal_effect_prediction_output_file=${expr_pred_output_dir}"marginal_eqtl_effect_prediction_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt.gz"
	ld_prune_r_sq_threshold="0.2"
	ld_prune_random_seed="1"
	ld_pruned_marginal_effect_prediction_output_file=${expr_pred_output_dir}"marginal_eqtl_effect_prediction_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_ld_pruned_r2_"${ld_prune_r_sq_threshold}"_seed_"${ld_prune_random_seed}".txt.gz"
	sbatch predict_marginal_eqtl_effects_from_borzoi_prior.sh $borzoi_effect_file $genotype_stem $genotype_sample_mapping_file $expr_file $borzoi_based_prior_output_stem $marginal_effect_prediction_output_file $ld_pruned_marginal_effect_prediction_output_file $ld_prune_r_sq_threshold $ld_prune_random_seed
fi





#####################
# Visualize expression predictions
if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_expression_prediction.R ${expr_pred_output_dir} $visualize_expression_pred
fi




#################
# Hybrid prediction stuff

target_tissue="Whole_Blood"
target_sample="GTEX-1LB8K-0005-SM-DIPED.1"
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"
max_training_sample_size="300"


	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${target_tissue}"_"${target_sample}"_"${anno_method}
if false; then
for training_sample_size in 25 50 75 100 125 150 175 200 225 250 300; do
	expr_prediction_output_file=${new_hybrid_expression_prediction_dir}"expr_pred_summary_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_"${training_sample_size}"_"${max_training_sample_size}".txt"
	sbatch hybrid_predict_cross_individual_expression.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno $training_sample_size $max_training_sample_size
	echo $training_sample_size
done
fi


if false; then
target_tissue="Whole_Blood"
target_sample="GTEX-1LB8K-0005-SM-DIPED.1"
ref_tissue="Muscle_Skeletal"
ref_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"
max_training_sample_size="300"


	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}


for training_sample_size in 25 50 75 100 125 150 175 200 225 250 300; do
	expr_prediction_output_file=${new_hybrid_expression_prediction_dir}"expr_pred_summary_"${target_tissue}"_"${target_sample}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_"${training_sample_size}"_"${max_training_sample_size}".txt"
	sbatch hybrid_predict_cross_individual_expression.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_prior_output_stem $expr_prediction_output_file $standardize_geno $training_sample_size $max_training_sample_size
done
fi



#####################
# Visualize expression predictions
if false; then
source ~/.bashrc
conda activate plink_env

Rscript visualize_expression_prediction.R ${expr_pred_output_dir} $visualize_expression_pred ${new_hybrid_expression_prediction_dir}
fi





#################
# Hybrid coloc thing
target_tissue="Liver"
target_sample="GTEX-11EQ9-0526-SM-5A5JZ.1"

anno_method="borzoi_effect_sizes"
standardize_geno="False"

trait_name="biochemistry_LDLdirect"
trait_ss_file=${gwas_sumstats_dir}${trait_name}"_hg38_liftover_sumstats.bgen.stats.gz"
if false; then
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	borzoi_annotation_file="NA"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

	training_sample_size="100"

	coloc_output_file=${hybrid_coloc_dir}"coloc_summary_"$trait_name"_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${standardize_geno}"_"${training_sample_size}
	sbatch hybrid_coloc.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $coloc_output_file $standardize_geno $training_sample_size $trait_ss_file


	training_sample_size="200"

	coloc_output_file=${hybrid_coloc_dir}"coloc_summary_"$trait_name"_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${standardize_geno}"_"${training_sample_size}
	sbatch hybrid_coloc.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $coloc_output_file $standardize_geno $training_sample_size $trait_ss_file
	training_sample_size="300"

	coloc_output_file=${hybrid_coloc_dir}"coloc_summary_"$trait_name"_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${standardize_geno}"_"${training_sample_size}
	sbatch hybrid_coloc.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $coloc_output_file $standardize_geno $training_sample_size $trait_ss_file

	training_sample_size="400"

	coloc_output_file=${hybrid_coloc_dir}"coloc_summary_"$trait_name"_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${standardize_geno}"_"${training_sample_size}
	sbatch hybrid_coloc.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $coloc_output_file $standardize_geno $training_sample_size $trait_ss_file
fi



#################
# TWAS
if false; then
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"

ref_tissue="Muscle_Skeletal"
ref_sample="GTEX-13QJ3-0726-SM-5SI68.1"

genotype_stem=${processed_genotype_data_dir}"gtex_v9_eqtl_chr"
borzoi_annotation_file="NA"
borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}

# trait_name<TAB>target_tissue<TAB>target_sample
while IFS=$'\t' read -r trait_name target_tissue target_sample; do

    echo "Submitting: ${trait_name} ${target_tissue} ${target_sample}"

    trait_ss_file=${gwas_sumstats_dir}${trait_name}"_hg38_liftover_sumstats.bgen.stats.gz"

    borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
    genotype_sample_mapping_file=${processed_genotype_data_dir}"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"
    expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

    twas_output_file=${twas_dir}"twas_summary_"${trait_name}"_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt"

    if [ ! -f "$trait_ss_file" ]; then
        echo "Missing trait sumstats: $trait_ss_file"
        continue
    fi

    if [ ! -f "$borzoi_effect_file" ]; then
        echo "Missing Borzoi effect file: $borzoi_effect_file"
        continue
    fi

    if [ ! -f "$genotype_sample_mapping_file" ]; then
        echo "Missing genotype sample mapping file: $genotype_sample_mapping_file"
        continue
    fi

    if [ ! -f "$expr_file" ]; then
        echo "Missing expression file: $expr_file"
        continue
    fi

    sbatch twas_run.sh \
        "$borzoi_effect_file" \
        "$borzoi_annotation_file" \
        "$genotype_stem" \
        "$genotype_sample_mapping_file" \
        "$expr_file" \
        "$distribution" \
        "$borzoi_based_prior_output_stem" \
        "$twas_output_file" \
        "$standardize_geno" \
        "$trait_ss_file"

done <<EOF
disease_AID_ALL	Spleen	GTEX-14PJ4-0526-SM-6871G.1
blood_MONOCYTE_COUNT	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
blood_MEAN_CORPUSCULAR_HEMOGLOBIN	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
blood_MEAN_PLATELET_VOL	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT	Cells_Cultured_fibroblasts	GTEX-139TS-0008-SM-62LDG.1
disease_ALLERGY_ECZEMA_DIAGNOSED	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
biochemistry_VitaminD	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
biochemistry_Cholesterol	Liver	GTEX-11EQ9-0526-SM-5A5JZ.1
repro_MENARCHE_AGE	Pituitary	GTEX-12WSC-3126-SM-5GCNB.1
bp_DIASTOLICadjMEDz	Artery_Aorta	GTEX-1JK1U-0426-SM-CYPSP.1
lung_FVCzSMOKE	Lung	GTEX-1399S-1726-SM-5L3DI.1
lung_FEV1FVCzSMOKE	Lung	GTEX-1399S-1726-SM-5L3DI.1
body_HEIGHTz	Cells_Cultured_fibroblasts	GTEX-139TS-0008-SM-62LDG.1
bmd_HEEL_TSCOREz	Cells_Cultured_fibroblasts	GTEX-139TS-0008-SM-62LDG.1
EOF
fi


# Non-gtex
if false; then
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"

# Fixed reference tissue for prior
ref_tissue="Muscle_Skeletal"
ref_sample="GTEX-13QJ3-0726-SM-5SI68.1"

borzoi_annotation_file="NA"
genotype_stem=${processed_genotype_data_dir}"gtex_v9_eqtl_chr"

borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}

# Columns:
# trait_name  non_gtex_tissue  non_gtex_sample  target_gtex_tissue  target_gtex_sample
while IFS=$'\t' read -r trait_name non_gtex_tissue non_gtex_sample target_gtex_tissue target_gtex_sample; do

    echo "Submitting: ${trait_name} | ${non_gtex_tissue} (${non_gtex_sample}) | GTEx proxy: ${target_gtex_tissue}"

    trait_ss_file=${gwas_sumstats_dir}${trait_name}"_hg38_liftover_sumstats.bgen.stats.gz"

    borzoi_effect_file=${borzoi_output_dir}${non_gtex_tissue}"_"${non_gtex_sample}"_borzoi_effects.txt.gz"
    genotype_sample_mapping_file=${processed_genotype_data_dir}"genotype_sample_mapping_to_"${target_gtex_tissue}"_expression_samples.txt"
    expr_file=${gtex_expr_dir}${target_gtex_tissue}".v10.residualized_expression_renormalized.bed"

    twas_output_file=${twas_dir}"twas_summary_"${trait_name}"_"${non_gtex_tissue}"_"${target_gtex_tissue}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt"

    sbatch twas_run.sh \
        "$borzoi_effect_file" \
        "$borzoi_annotation_file" \
        "$genotype_stem" \
        "$genotype_sample_mapping_file" \
        "$expr_file" \
        "$distribution" \
        "$borzoi_based_prior_output_stem" \
        "$twas_output_file" \
        "$standardize_geno" \
        "$trait_ss_file"


    borzoi_effect_file=${borzoi_output_dir}${target_gtex_tissue}"_"${target_gtex_sample}"_borzoi_effects.txt.gz"
    twas_output_file=${twas_dir}"twas_summary_"${trait_name}"_matched_proxy_for_"${non_gtex_tissue}"_"${target_gtex_tissue}"_"${ref_tissue}"_"${anno_method}"_"${distribution}"_"${standardize_geno}".txt"

    sbatch twas_run.sh \
        "$borzoi_effect_file" \
        "$borzoi_annotation_file" \
        "$genotype_stem" \
        "$genotype_sample_mapping_file" \
        "$expr_file" \
        "$distribution" \
        "$borzoi_based_prior_output_stem" \
        "$twas_output_file" \
        "$standardize_geno" \
        "$trait_ss_file"


done <<EOF
bmd_HEEL_TSCOREz	osteoblast	ENCFF953CRN	Muscle_Skeletal	GTEX-13QJ3-0726-SM-5SI68.1
bmd_HEEL_TSCOREz	osteocyte	ENCFF092QIW	Muscle_Skeletal	GTEX-13QJ3-0726-SM-5SI68.1
body_HEIGHTz	knee_articular_chondrocyte	ENCFF314QIG	Muscle_Skeletal	GTEX-13QJ3-0726-SM-5SI68.1
body_HEIGHTz	chondrocyte	ENCFF958WEW	Muscle_Skeletal	GTEX-13QJ3-0726-SM-5SI68.1
body_HEIGHTz	fetal_skeletal_muscle	ENCFF285OHV	Muscle_Skeletal	GTEX-13QJ3-0726-SM-5SI68.1
body_HEIGHTz	skeletal_muscle_satellite_cell	ENCFF804YRV	Muscle_Skeletal	GTEX-13QJ3-0726-SM-5SI68.1
disease_ASTHMA_DIAGNOSED	bronchial_epithelial_cell	ENCFF153YEN	Lung	GTEX-1399S-1726-SM-5L3DI.1
disease_ASTHMA_DIAGNOSED	bronchial_smooth_muscle_cell	ENCFF383NTW	Lung	GTEX-1399S-1726-SM-5L3DI.1
lung_FEV1FVCzSMOKE	alveolar_epithelial_cell	ENCFF734VZN	Lung	GTEX-1399S-1726-SM-5L3DI.1
lung_FVCzSMOKE	airway_epithelial_cell	ENCFF595CZA	Lung	GTEX-1399S-1726-SM-5L3DI.1
disease_RESPIRATORY_ENT	bronchus_fibroblast	ENCFF383RPD	Lung	GTEX-1399S-1726-SM-5L3DI.1
lung_FEV1FVCzSMOKE	fetal_lung	ENCFF892OBT	Lung	GTEX-1399S-1726-SM-5L3DI.1
disease_T2D	pancreatic_beta_cell	ENCFF995AUL	Pancreas	GTEX-11I78-0626-SM-5A5LZ.1
biochemistry_HbA1c	pancreatic_beta_cell	ENCFF995AUL	Pancreas	GTEX-11I78-0626-SM-5A5LZ.1
biochemistry_Glucose	endocrine_pancreas_progenitor_cell	ENCFF225CYB	Pancreas	GTEX-11I78-0626-SM-5A5LZ.1
disease_T2D	endocrine_pancreas_progenitor_cell	ENCFF225CYB	Pancreas	GTEX-11I78-0626-SM-5A5LZ.1
bp_SYSTOLICadjMEDz	aortic_smooth_muscle_cell	ENCFF281BWX	Artery_Aorta	GTEX-1JK1U-0426-SM-CYPSP.1
bp_DIASTOLICadjMEDz	aortic_smooth_muscle_cell	ENCFF281BWX	Artery_Aorta	GTEX-1JK1U-0426-SM-CYPSP.1
disease_HYPERTENSION_DIAGNOSED	aortic_smooth_muscle_cell	ENCFF281BWX	Artery_Aorta	GTEX-1JK1U-0426-SM-CYPSP.1
disease_CARDIOVASCULAR	coronary_artery_smooth_muscle_cell	ENCFF537AIY	Artery_Aorta	GTEX-1JK1U-0426-SM-CYPSP.1
disease_CARDIOVASCULAR	coronary_artery_endothelial_cell	ENCFF226UWU	Artery_Aorta	GTEX-1JK1U-0426-SM-CYPSP.1
bp_SYSTOLICadjMEDz	fetal_metanephros	ENCFF367PUX	Kidney_Cortex	GTEX-13112-2126-SM-5GCO4.1
disease_AID_ALL	regulatory_t_cell	ENCFF772RIQ	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
disease_PSORIASIS	activated_cd4_t_cell	ENCFF362PNM	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
disease_ALLERGY_ECZEMA_DIAGNOSED	keratinocyte	ENCFF703WGS	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
repro_MENARCHE_AGE	fetal_diencephalon	ENCFF355IUO	Pituitary	GTEX-12WSC-3126-SM-5GCNB.1
cov_EDU_COLLEGE	fetal_frontal_cortex	ENCFF217HQN	Brain_Cortex	GTEX-1H3O1-1726-SM-9WYSR.1
mental_NEUROTICISM	ganglionic_eminence_neurosphere	ENCFF674RUW	Brain_Cortex	GTEX-1H3O1-1726-SM-9WYSR.1
pigment_TANNING	adult_skin_melanocyte	ENCFF230ZSG	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
body_BALDING1	hair_follicle_dermal_papilla_cell	ENCFF798NEI	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
bp_SYSTOLICadjMEDz	glomerular_endothelial_cell	ENCFF863JIL	Kidney_Cortex	GTEX-13112-2126-SM-5GCO4.1
disease_HYPERTENSION_DIAGNOSED	mesangial_cell	ENCFF346QDJ	Kidney_Cortex	GTEX-13112-2126-SM-5GCO4.1
bp_DIASTOLICadjMEDz	proximal_tubule_epithelial_cell	ENCFF767MLU	Kidney_Cortex	GTEX-13112-2126-SM-5GCO4.1
repro_NumberChildrenEverBorn_Pooled	placental_epithelial_cell	ENCFF591NKN	Uterus	GTEX-13FTX-1026-SM-5J2O5.1
repro_NumberChildrenEverBorn_Pooled	placental_pericyte	ENCFF779FMR	Uterus	GTEX-13FTX-1026-SM-5J2O5.1
repro_NumberChildrenEverBorn_Pooled	villous_mesenchyme_fibroblast	ENCFF649CUP	Uterus	GTEX-13FTX-1026-SM-5J2O5.1
cov_EDU_COLLEGE	neural_progenitor_cell	ENCFF403DKN	Brain_Cortex	GTEX-1H3O1-1726-SM-9WYSR.1
mental_NEUROTICISM	cortex_neurosphere	ENCFF615IJV	Brain_Cortex	GTEX-1H3O1-1726-SM-9WYSR.1
other_MORNINGPERSON	fetal_diencephalon	ENCFF355IUO	Pituitary	GTEX-12WSC-3126-SM-5GCNB.1
body_HEIGHTz	neural_crest_cell	ENCFF024MEB	Brain_Cortex	GTEX-1H3O1-1726-SM-9WYSR.1
disease_ALLERGY_ECZEMA_DIAGNOSED	dermal_fibroblast	ENCFF343FWG	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
body_BALDING1	hair_follicular_keratinocyte	ENCFF544MJM	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
disease_CARDIOVASCULAR	aortic_adventitial_fibroblast	ENCFF256ZZS	Artery_Aorta	GTEX-1JK1U-0426-SM-CYPSP.1
disease_CARDIOVASCULAR	cardiac_ventricle_fibroblast	ENCFF287LRZ	Heart_Left_Ventricle	GTEX-18465-0926-SM-731AY.1
disease_AID_ALL	activated_t_cell_il2_cd3_cd28	ENCFF122XJD	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
disease_PSORIASIS	activated_cd4_memory_t_cell	ENCFF330EUN	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
EOF

fi



# TWAS Downsample
if false; then
anno_method="borzoi_effect_sizes"
distribution="ldscore_grid_squared"
standardize_geno="False"

ref_tissue="Muscle_Skeletal"
ref_sample="GTEX-13QJ3-0726-SM-5SI68.1"

genotype_stem=${processed_genotype_data_dir}"gtex_v9_eqtl_chr"
borzoi_annotation_file="NA"
borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}

for training_sample_size in 25 50 100 150 200 250; do

# trait_name<TAB>target_tissue<TAB>target_sample
while IFS=$'\t' read -r trait_name target_tissue target_sample; do

    echo "Submitting: ${trait_name} ${target_tissue} ${target_sample}"

    trait_ss_file=${gwas_sumstats_dir}${trait_name}"_hg38_liftover_sumstats.bgen.stats.gz"

    borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
    genotype_sample_mapping_file=${processed_genotype_data_dir}"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"
    expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"

    twas_output_file=${twas_dir}"twas_summary_"${trait_name}"_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_downsample_"${training_sample_size}".txt"

    if [ ! -f "$trait_ss_file" ]; then
        echo "Missing trait sumstats: $trait_ss_file"
        continue
    fi

    if [ ! -f "$borzoi_effect_file" ]; then
        echo "Missing Borzoi effect file: $borzoi_effect_file"
        continue
    fi

    if [ ! -f "$genotype_sample_mapping_file" ]; then
        echo "Missing genotype sample mapping file: $genotype_sample_mapping_file"
        continue
    fi

    if [ ! -f "$expr_file" ]; then
        echo "Missing expression file: $expr_file"
        continue
    fi

    sbatch twas_run_downsample.sh \
        "$borzoi_effect_file" \
        "$borzoi_annotation_file" \
        "$genotype_stem" \
        "$genotype_sample_mapping_file" \
        "$expr_file" \
        "$distribution" \
        "$borzoi_based_prior_output_stem" \
        "$twas_output_file" \
        "$standardize_geno" \
        "$trait_ss_file" \
        "$training_sample_size"

done <<EOF
disease_AID_ALL	Spleen	GTEX-14PJ4-0526-SM-6871G.1
blood_MONOCYTE_COUNT	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
blood_MEAN_CORPUSCULAR_HEMOGLOBIN	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
blood_MEAN_PLATELET_VOL	Whole_Blood	GTEX-1LB8K-0005-SM-DIPED.1
blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT	Cells_Cultured_fibroblasts	GTEX-139TS-0008-SM-62LDG.1
disease_ALLERGY_ECZEMA_DIAGNOSED	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
biochemistry_VitaminD	Skin_Sun_Exposed_Lower_leg	GTEX-13U4I-0126-SM-5LU38.1
biochemistry_Cholesterol	Liver	GTEX-11EQ9-0526-SM-5A5JZ.1
repro_MENARCHE_AGE	Pituitary	GTEX-12WSC-3126-SM-5GCNB.1
bp_DIASTOLICadjMEDz	Artery_Aorta	GTEX-1JK1U-0426-SM-CYPSP.1
lung_FVCzSMOKE	Lung	GTEX-1399S-1726-SM-5L3DI.1
lung_FEV1FVCzSMOKE	Lung	GTEX-1399S-1726-SM-5L3DI.1
body_HEIGHTz	Cells_Cultured_fibroblasts	GTEX-139TS-0008-SM-62LDG.1
bmd_HEEL_TSCOREz	Cells_Cultured_fibroblasts	GTEX-139TS-0008-SM-62LDG.1
EOF

done
fi




if false; then
source ~/.bashrc
conda activate plink_env

trait_name="bmd_HEEL_TSCOREz"
non_gtex_tissue="osteoblast"
non_gtex_sample="ENCFF953CRN"
ref_tissue="Muscle_Skeletal"
ref_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_effect_sizes"
gene_id="ENSG00000152779"

borzoi_effect_file=${borzoi_output_dir}${non_gtex_tissue}"_"${non_gtex_sample}"_borzoi_effects.txt.gz"
genotype_stem=${processed_genotype_data_dir}"gtex_v9_eqtl_chr"
trait_ss_file=${gwas_sumstats_dir}${trait_name}"_hg38_liftover_sumstats.bgen.stats.gz"
borzoi_based_prior_output_stem=${borzoi_based_prior_output_dir}"borzoi_based_prior_ldscore_grid_squared_based_"${ref_tissue}"_"${ref_sample}"_"${anno_method}
specific_example_output_file=${twas_dir}"specific_example_data_"${trait_name}"_"${non_gtex_tissue}"_"${gene_id}".txt"
specific_example_plot_output_file=${twas_dir}"specific_example_locus_plot_"${trait_name}"_"${non_gtex_tissue}"_"${gene_id}".pdf"

python extract_specific_example_data.py \
	"$borzoi_effect_file" \
	"$genotype_stem" \
	"$trait_ss_file" \
	"$borzoi_based_prior_output_stem" \
	"$gene_id" \
	"$specific_example_output_file"

python plot_specific_example_locus.py \
	"$specific_example_output_file" \
	"$specific_example_plot_output_file"

echo $specific_example_output_file
echo $specific_example_plot_output_file
fi


if false; then
source ~/.bashrc
conda activate plink_env
python organize_twas_results.py ${twas_dir}

Rscript visualize_twas_results.R ${twas_dir}
fi





################################################
# Kinda a seperate track on hybrid prediction of expression (compared to elastic net)
################################################

#################
# Borzoi gaussian prior

target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="gaussian"
standardize_geno="True"
n_folds="5"
down_sampling_fraction="1.0"


borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"
if false; then
for fold_number in $(seq 1 ${n_folds}); do
	borzoi_based_hybrid_expr_pred_output_stem=${hybrid_expression_prediction_dir}"borzoi_based_hybrid_expr_pred_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_"${fold_number}"_"${n_folds}"_"${down_sampling_fraction}
	sbatch run_borzoi_based_hybrid_expression_prediction.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_hybrid_expr_pred_output_stem $standardize_geno $fold_number $n_folds $down_sampling_fraction
done
fi


#################
# Elastic net
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="elastic_net"
standardize_geno="True"
n_folds="5"
down_sampling_fraction="1.0"


borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"
if false; then
for fold_number in $(seq 1 ${n_folds}); do
	borzoi_based_hybrid_expr_pred_output_stem=${hybrid_expression_prediction_dir}"borzoi_based_hybrid_expr_pred_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_"${fold_number}"_"${n_folds}"_"${down_sampling_fraction}
	sbatch run_borzoi_based_hybrid_expression_prediction.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_hybrid_expr_pred_output_stem $standardize_geno $fold_number $n_folds $down_sampling_fraction
done
fi

#################
# Genome-wide  no borzoi ridge
target_tissue="Muscle_Skeletal"
target_sample="GTEX-13QJ3-0726-SM-5SI68.1"
anno_method="borzoi_magnitude_bins"
distribution="no_borzoi_genome_wide_ridge"
standardize_geno="True"
n_folds="5"
down_sampling_fraction="1.0"


borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_"${anno_method}"_annotations.txt.gz"
genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

expr_file=${gtex_expr_dir}${target_tissue}".v10.residualized_expression_renormalized.bed"
if false; then
for fold_number in $(seq 1 ${n_folds}); do
	borzoi_based_hybrid_expr_pred_output_stem=${hybrid_expression_prediction_dir}"borzoi_based_hybrid_expr_pred_results_"${target_tissue}"_"${target_sample}"_"${anno_method}"_"${distribution}"_"${standardize_geno}"_"${fold_number}"_"${n_folds}"_"${down_sampling_fraction}
	sbatch run_borzoi_based_hybrid_expression_prediction.sh $borzoi_effect_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file $expr_file $distribution $borzoi_based_hybrid_expr_pred_output_stem $standardize_geno $fold_number $n_folds $down_sampling_fraction
done
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
