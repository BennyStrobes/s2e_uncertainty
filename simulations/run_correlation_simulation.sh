#!/bin/bash
#SBATCH -t 0-20:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                          # Partition to run in
#SBATCH --mem=20GB  



simulation_iter="${1}"
gene_ld_summary_file="${2}"
causal_effect_dir="${3}"
est_eqtl_effect_size_dir="${4}"
est_borzoi_effect_size_dir="${5}"
onek_genomes_plink_filestem="${6}"
inf_output_dir="${7}"
sldmc_code_dir="${8}"

echo "Simulation "${simulation_iter}
source ~/.bashrc
conda activate plink_env


date
####################################################
# Part 1: Simulate causal variant-gene effect sizes
####################################################
echo "PART 1"
causal_variant_gene_effect_size_file=${causal_effect_dir}"sim"${simulation_iter}"_sim_causal_variant_gene_effects.txt.gz"
if false; then
python simulate_causal_variant_gene_effect_size.py $simulation_iter $gene_ld_summary_file $causal_variant_gene_effect_size_file
fi


####################################################
# Part 2: Simulate estimated borzoi (standardized effect) sizes
####################################################
echo "PART 2"
n_anno="6"
est_borzoi_effect_size_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_est_borzoi_effects_"${n_anno}"_anno.txt.gz"
sim_variant_gene_annotation_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_variant_gene_annotations_"${n_anno}"_annotations.txt.gz"
# SLDMC derives the category file name from the annotation file name, so these two must stay paired
sldmc_variant_gene_annotation_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_variant_gene_annotations_"${n_anno}"_annotations_sldmc.txt.gz"
sldmc_annotation_category_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_variant_gene_annotations_"${n_anno}"_annotations_sldmc_categories.txt"
if false; then
python simulate_est_borzoi_effects_for_correlation_experiment.py $causal_variant_gene_effect_size_file $est_borzoi_effect_size_file ${simulation_iter} $sim_variant_gene_annotation_file $n_anno $sldmc_variant_gene_annotation_file $sldmc_annotation_category_file
fi

####################################################
# Part 2.5: Generate true simulated calibration effect sizes + correlation
####################################################
echo "PART 2.5"
simulation_parameter_summary_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_variant_gene_annotations_"${n_anno}"_true_sim_effect_summary.txt"
if false; then
python calculate_true_simulated_calibration_effect_sizes_and_correlation.py $est_borzoi_effect_size_file $sim_variant_gene_annotation_file $causal_variant_gene_effect_size_file $simulation_parameter_summary_file
fi

####################################################
# Part 3: Simulate estimated eqtl effect sizes
####################################################
echo "PART 3"
eqtl_sample_size="489"
est_eqtl_effect_size_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_est_eqtl_effects.txt.gz"
ind_expr_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_individual_expression.txt.gz"
susie_fine_mapping_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_susie_fine_mapping.txt.gz"
genotype_sample_mapping_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_genotype_sample_mapping.txt"
source ~/.bashrc
conda activate susie
python simulate_eqtl_analysis.py $causal_variant_gene_effect_size_file $est_eqtl_effect_size_file $gene_ld_summary_file $onek_genomes_plink_filestem $eqtl_sample_size $simulation_iter $ind_expr_file $susie_fine_mapping_file $genotype_sample_mapping_file




####################################################
# Part 4: Run LD corr inference
####################################################
if false; then
echo "PART 4"
# Updated code
source ~/.bashrc
conda activate sldmc
ld_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_ld_corr_results"
python ${sldmc_code_dir}sldmc.py \
    --est-borzoi-effect-size-file $est_borzoi_effect_size_file \
    --est-eqtl-effect-size-file $est_eqtl_effect_size_file \
    --sim-variant-gene-annotation-file $sldmc_variant_gene_annotation_file \
    --genotype-plink-filestem $onek_genomes_plink_filestem \
    --genotype-sample-mapping-file $genotype_sample_mapping_file \
    --ld-corr-output-stem $ld_corr_output_stem 
fi


if false; then
# OLD Version of this
source ~/.bashrc
conda activate plink_env
echo "PART 4"
ld_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_ld_corr_results"
python run_ld_corr.py $est_borzoi_effect_size_file $est_eqtl_effect_size_file $sim_variant_gene_annotation_file $onek_genomes_plink_filestem $ld_corr_output_stem
fi


####################################################
# Part 5: Run correlations based on only fine-mapped snps
####################################################
if false; then
source ~/.bashrc
conda activate plink_env
echo "PART 5"
fm_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_fm_corr_results"
python run_fine_map_corr.py $est_borzoi_effect_size_file $susie_fine_mapping_file $sim_variant_gene_annotation_file $onek_genomes_plink_filestem $fm_corr_output_stem
fi




date
