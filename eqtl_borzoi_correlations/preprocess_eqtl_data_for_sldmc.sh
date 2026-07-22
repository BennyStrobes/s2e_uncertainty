#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 





gtex_v10_pc_genes_gtf="${1}"
eqtl_sumstats_dir="${2}"
tissue_name="${3}"
eqtl_ss_output_dir="${4}"



source ~/.bashrc
conda activate plink_env


echo $tissue_name

python preprocess_eqtl_data_for_sldmc.py $gtex_v10_pc_genes_gtf $eqtl_sumstats_dir $tissue_name $eqtl_ss_output_dir