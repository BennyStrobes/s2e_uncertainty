#!/bin/bash
#SBATCH --gpus 1                             # Request one core
#SBATCH -t 0-50:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-gpu                           # Partition to run in
#SBATCH --mem=30GB                         # Memory total in MiB (for all cores)
#SBATCH --cpus-per-task=8                        # Number of CPUs per task

source ~/.bashrc
conda activate borzoi


. /programs/local/miniforge/etc/profile.d/conda.sh; conda activate borzoi; echo $HOSTNAME; hound_train.py      -o /lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_testing/tutorial/model_training/micro_models/f0c0/train     /lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_testing/tutorial/model_training/micro_models/f0c0/params.json /lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/borzoi_testing/tutorial/model_training/micro_models/f0c0/data0