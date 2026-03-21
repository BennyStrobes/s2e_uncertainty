#!/bin/bash
#SBATCH -t 0-3:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                          # Partition to run in
#SBATCH --mem=10GB  

######################
# Input data
######################
# File stem to 1 thousand genomes plink files
onek_genomes_plink_filestem="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/1000G.EUR.hg38."

# Gene annotation file
gene_annotation_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"




######################
# Output data
######################
# Output root
output_root="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/simulations/"

# LD directory
gene_info_dir=${output_root}"gene_info/"

# LD directory
ld_dir=${output_root}"LD/"

# simulated causal eqtl effects directory
causal_effect_dir=${output_root}"sim_causal_eqtl_effects/"

# simulated estimated eqtl effect size directory
est_eqtl_effect_size_dir=${output_root}"sim_est_eqtl_effects/"

# simulated estimated borzoi effect size directory
est_borzoi_effect_size_dir=${output_root}"sim_est_borzoi_effects/"

# directory containing ldsc-like estimates of sampling variance
ldsc_like_samp_var_dir=${output_root}"ldsc_like_samp_var/"


# directory containing ashr-like estimates of sampling variance
ashr_like_samp_var_dir=${output_root}"ashr_like_samp_var/"

# directory containing visualizations fo results
visualization_dir=${output_root}"visualization_dir/"


######################
# Run analysis
######################

# First generate gene info including list of genes and variant-to-gene links
if false; then
source ~/.bashrc
conda activate plink_env
python generate_simulation_gene_info.py $gene_annotation_file $onek_genomes_plink_filestem $gene_info_dir
fi


# Second generate LD data
gene_summary_file=$gene_info_dir"gene_summary.txt"
if false; then
source ~/.bashrc
conda activate plink_env
python generate_ld_scores_and_ld_mean_for_simulation.py $onek_genomes_plink_filestem $gene_summary_file $ld_dir
fi

gene_ld_summary_file=${ld_dir}"gene_summary_w_ld_full2.txt"



######################
# Run simulations
######################
if false; then
for simulation_iter in {1..50}
do
	sbatch run_simulation.sh $simulation_iter $gene_ld_summary_file $causal_effect_dir $est_eqtl_effect_size_dir $est_borzoi_effect_size_dir $onek_genomes_plink_filestem $ldsc_like_samp_var_dir $ashr_like_samp_var_dir
done
fi

if false; then
simulation_iter="1"
sbatch run_simulation.sh $simulation_iter $gene_ld_summary_file $causal_effect_dir $est_eqtl_effect_size_dir $est_borzoi_effect_size_dir $onek_genomes_plink_filestem $ldsc_like_samp_var_dir $ashr_like_samp_var_dir

simulation_iter="2"
sbatch run_simulation.sh $simulation_iter $gene_ld_summary_file $causal_effect_dir $est_eqtl_effect_size_dir $est_borzoi_effect_size_dir $onek_genomes_plink_filestem $ldsc_like_samp_var_dir $ashr_like_samp_var_dir

simulation_iter="3"
sbatch run_simulation.sh $simulation_iter $gene_ld_summary_file $causal_effect_dir $est_eqtl_effect_size_dir $est_borzoi_effect_size_dir $onek_genomes_plink_filestem $ldsc_like_samp_var_dir $ashr_like_samp_var_dir

simulation_iter="4"
sbatch run_simulation.sh $simulation_iter $gene_ld_summary_file $causal_effect_dir $est_eqtl_effect_size_dir $est_borzoi_effect_size_dir $onek_genomes_plink_filestem $ldsc_like_samp_var_dir $ashr_like_samp_var_dir

simulation_iter="5"
sbatch run_simulation.sh $simulation_iter $gene_ld_summary_file $causal_effect_dir $est_eqtl_effect_size_dir $est_borzoi_effect_size_dir $onek_genomes_plink_filestem $ldsc_like_samp_var_dir $ashr_like_samp_var_dir
fi



simulation_iter="6"
if false; then
sbatch run_simulation.sh $simulation_iter $gene_ld_summary_file $causal_effect_dir $est_eqtl_effect_size_dir $est_borzoi_effect_size_dir $onek_genomes_plink_filestem $ldsc_like_samp_var_dir $ashr_like_samp_var_dir
fi






if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_simulation_results.R $ldsc_like_samp_var_dir $visualization_dir
fi



