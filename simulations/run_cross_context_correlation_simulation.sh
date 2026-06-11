#!/bin/bash
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                          # Partition to run in
#SBATCH --mem=30GB  



simulation_iter="${1}"
gene_ld_summary_file="${2}"
causal_effect_dir="${3}"
est_eqtl_effect_size_dir="${4}"
est_borzoi_effect_size_dir="${5}"
onek_genomes_plink_filestem="${6}"
inf_output_dir="${7}"
rho_delta="${8}"

echo "Simulation "${simulation_iter}
source ~/.bashrc
conda activate plink_env


date
####################################################
# Part 1: Simulate causal variant-gene effect sizes
####################################################
echo "PART 1"
t1_causal_variant_gene_effect_size_file=${causal_effect_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_t1_causal_variant_gene_effects.txt.gz"
t2_causal_variant_gene_effect_size_file=${causal_effect_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_t2_causal_variant_gene_effects.txt.gz"
if false; then
python simulate_xc_causal_variant_gene_effect_size.py $simulation_iter $gene_ld_summary_file $t1_causal_variant_gene_effect_size_file $t2_causal_variant_gene_effect_size_file
fi

####################################################
# Part 2: Simulate estimated borzoi (standardized effect) sizes
####################################################
echo "PART 2"
t1_est_borzoi_effect_size_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_t1_est_borzoi_effects.txt.gz"
t2_est_borzoi_effect_size_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_t2_est_borzoi_effects.txt.gz"
if false; then
python simulate_xc_est_borzoi_effects_for_correlation_experiment.py $t1_causal_variant_gene_effect_size_file $t2_causal_variant_gene_effect_size_file $t1_est_borzoi_effect_size_file $t2_est_borzoi_effect_size_file ${simulation_iter} $rho_delta
fi

####################################################
# Part 2.5: Generate true simulated calibration effect sizes + correlation
####################################################
echo "PART 2.5"
simulation_parameter_summary_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_variant_gene_true_sim_effect_summary.txt"
if false; then
python calculate_xc_true_simulated_calibration_effect_sizes_and_correlation.py $t1_est_borzoi_effect_size_file $t2_est_borzoi_effect_size_file $t1_causal_variant_gene_effect_size_file $t2_causal_variant_gene_effect_size_file $simulation_parameter_summary_file
fi



####################################################
# Part 3: Simulate estimated eqtl effect sizes
####################################################
echo "PART 3: tissue 1"
t1_eqtl_sample_size="460"
t1_est_eqtl_effect_size_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_eqtl_ss_"${t1_eqtl_sample_size}"_t1_est_eqtl_effects.txt.gz"
t1_ind_expr_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_eqtl_ss_"${t1_eqtl_sample_size}"_t1_individual_expression.txt.gz"
t1_susie_fine_mapping_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_eqtl_ss_"${t1_eqtl_sample_size}"_t1_susie_fine_mapping.txt.gz"
if false; then
source ~/.bashrc
conda activate susie
python simulate_eqtl_analysis.py $t1_causal_variant_gene_effect_size_file $t1_est_eqtl_effect_size_file $gene_ld_summary_file $onek_genomes_plink_filestem $t1_eqtl_sample_size $simulation_iter $t1_ind_expr_file $t1_susie_fine_mapping_file
fi


echo "PART 3: tissue 2"
t2_eqtl_sample_size="390"
t2_simulation_iter=$((simulation_iter + 100000))
t2_est_eqtl_effect_size_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_eqtl_ss_"${t2_eqtl_sample_size}"_t2_est_eqtl_effects.txt.gz"
t2_ind_expr_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_eqtl_ss_"${t2_eqtl_sample_size}"_t2_individual_expression.txt.gz"
t2_susie_fine_mapping_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_sim_eqtl_ss_"${t2_eqtl_sample_size}"_t2_susie_fine_mapping.txt.gz"
if false; then
source ~/.bashrc
conda activate susie
python simulate_eqtl_analysis.py $t2_causal_variant_gene_effect_size_file $t2_est_eqtl_effect_size_file $gene_ld_summary_file $onek_genomes_plink_filestem $t2_eqtl_sample_size $t2_simulation_iter $t2_ind_expr_file $t2_susie_fine_mapping_file
fi


####################################################
# Part 4: Run cross-context LD corr inference
####################################################
source ~/.bashrc
conda activate plink_env
echo "PART 4"
ld_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_t1_eqtl_ss_"${t1_eqtl_sample_size}"_t2_eqtl_ss_"${t2_eqtl_sample_size}"_xc_ld_corr_results"
if false; then
python run_cross_tissue_ld_corr.py $t1_est_borzoi_effect_size_file $t2_est_borzoi_effect_size_file $t1_est_eqtl_effect_size_file $t2_est_eqtl_effect_size_file $onek_genomes_plink_filestem $t1_eqtl_sample_size $t2_eqtl_sample_size $ld_corr_output_stem
fi

####################################################
# Part 5: Run correlations based on only fine-mapped snps
####################################################
source ~/.bashrc
conda activate plink_env
echo "PART 5"
fm_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_rho_"${rho_delta}"_t1_eqtl_ss_"${t1_eqtl_sample_size}"_t2_eqtl_ss_"${t2_eqtl_sample_size}"_xc_fm_corr_results"
python run_cross_tissue_fine_map_corr.py $t1_est_borzoi_effect_size_file $t2_est_borzoi_effect_size_file $t1_susie_fine_mapping_file $t2_susie_fine_mapping_file $t1_causal_variant_gene_effect_size_file $t2_causal_variant_gene_effect_size_file $fm_corr_output_stem
