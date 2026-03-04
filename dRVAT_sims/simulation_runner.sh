#!/bin/bash
#SBATCH -t 0-00:40                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 


sim_number="${1}"
output_stem="${2}"






source ~/.bashrc
conda activate borzoi


sample_size="50000"


n_detected="0"
python dRVAT_simulation.py $sim_number $sample_size ${n_detected} $output_stem"_"${sample_size}"_"${n_detected}


if false; then
for n_detected in {1..5}
do
    python dRVAT_simulation.py $sim_number $sample_size ${n_detected} $output_stem"_"${sample_size}"_"${n_detected}
done
fi


if false; then
sample_size="100000"
for n_detected in {1..5}
do
    python dRVAT_simulation.py $sim_number $sample_size ${n_detected} $output_stem"_"${sample_size}"_"${n_detected}
done
fi