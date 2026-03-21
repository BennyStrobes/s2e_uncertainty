#!/bin/bash
#SBATCH -t 0-30:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                          # Partition to run in
#SBATCH --mem=20GB  


simulation_iter="${1}"
gene_ld_summary_file="${2}"
causal_effect_dir="${3}"
est_eqtl_effect_size_dir="${4}"
est_borzoi_effect_size_dir="${5}"
onek_genomes_plink_filestem="${6}"
ldsc_like_samp_var_dir="${7}"
ashr_like_samp_var_dir="${8}"

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
borzoi_error_distribution="simple_mog"
if false; then
for sampling_variance in 0.002 0.004 0.006 0.008
do
	est_borzoi_effect_size_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_"${borzoi_error_distribution}"_"${sampling_variance}"_est_borzoi_effects.txt.gz"
	sim_borzoi_error_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_"${borzoi_error_distribution}"_"${sampling_variance}"_sim_borzoi_sampling_error.txt.gz"
	python simulate_est_borzoi_effects.py $causal_variant_gene_effect_size_file $est_borzoi_effect_size_file $borzoi_error_distribution ${simulation_iter} ${sampling_variance} ${sim_borzoi_error_file}
done
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
# Part 4: Infer sampling variance using ldsc-like approach
####################################################
echo "PART 4"
eqtl_sample_size="489"
borzoi_error_distribution="simple_mog"
if false; then
for sampling_variance in 0.002 0.004 0.006 0.008
do

	est_borzoi_effect_size_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_"${borzoi_error_distribution}"_"${sampling_variance}"_est_borzoi_effects.txt.gz"
	est_eqtl_effect_size_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_est_eqtl_effects.txt.gz"

	ldsc_results_file=${ldsc_like_samp_var_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${borzoi_error_distribution}"_"${sampling_variance}"_ldsc_like_est_of_sampling_variance.txt"

	python infer_borzoi_sampling_variance_with_ldsc_like_approach.py ${est_borzoi_effect_size_file} ${est_eqtl_effect_size_file} $gene_ld_summary_file $onek_genomes_plink_filestem $eqtl_sample_size $ldsc_results_file
done
fi



####################################################
# Part 5: Infer sampling variance distribution using ashr style approach
####################################################
echo "PART 5"
eqtl_sample_size="489"
borzoi_error_distribution="simple_mog"
sampling_variance="0.002"

	est_borzoi_effect_size_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_"${borzoi_error_distribution}"_"${sampling_variance}"_est_borzoi_effects.txt.gz"
	est_eqtl_effect_size_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_est_eqtl_effects.txt.gz"

	ashr_results_file=${ashr_like_samp_var_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${borzoi_error_distribution}"_"${sampling_variance}"_ashr_like_est_of_sampling_variance.txt"
	if false; then
	python infer_borzoi_sampling_distribution_with_ashr_like_approach.py ${est_borzoi_effect_size_file} ${est_eqtl_effect_size_file} $gene_ld_summary_file $onek_genomes_plink_filestem $eqtl_sample_size $ashr_results_file
	fi


	causal_ashr_results_file=${ashr_like_samp_var_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${borzoi_error_distribution}"_"${sampling_variance}"_causal_ashr_like_est_of_sampling_variance.txt"
	if false; then
	python infer_borzoi_sampling_distribution_with_causal_ashr_like_approach.py ${est_borzoi_effect_size_file} ${est_eqtl_effect_size_file} $gene_ld_summary_file $onek_genomes_plink_filestem $eqtl_sample_size $causal_ashr_results_file
	fi

	causal_ashr_individual_level_results_file=${ashr_like_samp_var_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${borzoi_error_distribution}"_"${sampling_variance}"_causal_ashr_like_individual_likelihood_est_of_sampling_variance.txt"
	python infer_borzoi_sampling_distribution_with_causal_ashr_individual_like_approach.py ${est_borzoi_effect_size_file} ${ind_expr_file} $gene_ld_summary_file $onek_genomes_plink_filestem $eqtl_sample_size $causal_ashr_individual_level_results_file


	sim_borzoi_error_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_"${borzoi_error_distribution}"_"${sampling_variance}"_sim_borzoi_sampling_error.txt.gz"
	output_stem=${ashr_like_samp_var_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${borzoi_error_distribution}"_"${sampling_variance}"_causal_ashr_individual"
	python visualize_sampling_cdf.py $causal_ashr_individual_level_results_file $sim_borzoi_error_file $output_stem

date


