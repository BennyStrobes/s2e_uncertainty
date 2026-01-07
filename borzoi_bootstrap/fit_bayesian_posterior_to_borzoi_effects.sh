#!/bin/bash
#SBATCH -t 0-24:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                           # Partition to run in
#SBATCH --mem=10GB  



organized_borzoi_bs_effecs_file="${1}"
bayes_fit_output_root="${2}"

source ~/.bashrc
conda activate borzoi


python fit_bayesian_posterior_to_borzoi_effects.py $organized_borzoi_bs_effecs_file $bayes_fit_output_root