#!/bin/bash
#SBATCH -t 0-5:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 




borzoi_results_dir="${1}"
plink2_genotype_data_dir="${2}"
eqtl_sumstats_dir="${3}"
borzoi_output_dir="${4}"
LD_output_dir="${5}"
eqtl_ss_output_dir="${6}"
chrom_num="${7}"
gtex_v10_pc_genes_gtf="${8}"

source ~/.bashrc
conda activate plink_env


python integrate_eqtl_borzoi_ld_data.py ${borzoi_results_dir} ${plink2_genotype_data_dir} ${eqtl_sumstats_dir} ${borzoi_output_dir} ${LD_output_dir} ${eqtl_ss_output_dir} $chrom_num $gtex_v10_pc_genes_gtf

