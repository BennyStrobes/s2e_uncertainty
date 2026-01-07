#!/bin/bash
#SBATCH -t 0-15:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB 




sardinia_raw_sumstat_file="${1}"
fm_eqtl_sumstats_file="${2}"
pip_threshold="${3}"

source ~/.bashrc
conda activate borzoi

python fine_map_eqtl_sumstats.py $sardinia_raw_sumstat_file $fm_eqtl_sumstats_file $pip_threshold
