#######################
# Input data
#######################


# Directory containing genotype data
processed_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/plink_processed_genotype/"

# Directory containing eQTL summary statistics
eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/eqtl_results/"


gwas_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/sumstats/UKBB_all_snps_sumstats/data/"


#######################
# Output directors
#######################
# root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/eqtl_borzoi_correlations/twas/"

# Twas simulation dir
twas_simulation_dir=${output_root}"twas_simulations/"






#####################
# Code
#####################

# Twas simulation
sim_fsr="0.2"
alpha_var="5e-4"
gwas_ss="300000"

twas_simulation_output_stem=${twas_simulation_dir}"twas_sim_"${sim_fsr}"_"$alpha_var"_"$gwas_ss"_results"
if false; then
sh run_twas_simulation.sh $processed_genotype_data_dir $twas_simulation_output_stem $sim_fsr $alpha_var $gwas_ss
fi














