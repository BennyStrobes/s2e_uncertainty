#!/bin/bash
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB 


$processed_genotype_data_dir $twas_simulation_output_stem $sim_fsr $alpha_var $gwas_ss


processed_genotype_data_dir="${1}"
twas_simulation_output_stem="${2}"
sim_fsr="${3}"
alpha_var="${4}"
gwas_ss="${5}"


source ~/.bashrc
conda activate plink_env

echo $twas_simulation_output_stem


python run_noisy_twas_simulation.py $processed_genotype_data_dir $twas_simulation_output_stem $sim_fsr $alpha_var $gwas_ss