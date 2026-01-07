#!/bin/bash
#SBATCH -t 0-15:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB 




fm_eqtl_sumstats_file="${1}"
fm_and_af_eqtl_sumstats_file="${2}"
hg19_1000G_genotype_dir="${3}"

source ~/.bashrc
conda activate borzoi

python add_allele_frequency_info_to_fm_file.py $fm_eqtl_sumstats_file $fm_and_af_eqtl_sumstats_file $hg19_1000G_genotype_dir
