


#################
# Input data
##################


borzoi_eqtl_effects_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_tissue_bootstrap/borzoi_pred_eqtl_effects_old/cross_bootstrap_PIP_0.9_borzoi_pred_eqtl_effects_cross_tissue.txt"

###################
# Output directories
#####################
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/dRVAT_sims/"

# Raw simulation results
sim_results_dir=${output_root}"sim_results/"

# organized simulation results
organized_sim_results_dir=${output_root}"organized_sim_results/"

# organized simulation results
visualization_dir=${output_root}"visualization/"




#########################
# Run simulation
#########################
if false; then
for sim_number in {1..100}
do
	sbatch simulation_runner.sh $sim_number $sim_results_dir"sim_"${sim_number}
done
fi



source ~/.bashrc
conda activate borzoi
python organize_simulation_results.py $sim_results_dir $organized_sim_results_dir


if false; then
source ~/.bashrc
conda activate plink_env
fi
if false; then
Rscript visualize_dRVAT_results.R $organized_sim_results_dir $visualization_dir $borzoi_eqtl_effects_file
fi

