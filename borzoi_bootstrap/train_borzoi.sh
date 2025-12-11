#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-95:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu-pe                           # Partition to run in
#SBATCH --mem=10GB                         # Memory total in MiB (for all cores)
#SBATCH --cpus-per-task=4                        # Number of CPUs per task




source ~/.bashrc
conda activate borzoi

submission_command="${1}"

python $submission_command

