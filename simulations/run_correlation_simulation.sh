#!/bin/bash
#SBATCH -t 0-3:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                          # Partition to run in
#SBATCH --mem=20GB  



simulation_iter="${1}"
gene_ld_summary_file="${2}"
causal_effect_dir="${3}"
est_eqtl_effect_size_dir="${4}"
est_borzoi_effect_size_dir="${5}"
onek_genomes_plink_filestem="${6}"
inf_output_dir="${7}"

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
if false; then
python simulate_est_borzoi_effects_for_correlation_experiment.py $causal_variant_gene_effect_size_file $est_borzoi_effect_size_file ${simulation_iter} $sim_variant_gene_annotation_file $n_anno
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
if false; then
python simulate_eqtl_analysis.py $causal_variant_gene_effect_size_file $est_eqtl_effect_size_file $gene_ld_summary_file $onek_genomes_plink_filestem $eqtl_sample_size $simulation_iter $ind_expr_file
fi



####################################################
# Part 4: Simulate estimated eqtl effect sizes
####################################################
echo "PART 4"
ld_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_ld_corr_results"
python run_ld_corr.py $est_borzoi_effect_size_file $est_eqtl_effect_size_file $sim_variant_gene_annotation_file $onek_genomes_plink_filestem $ld_corr_output_stem












