#!/bin/bash --login
#SBATCH --job-name=LTIFM
#SBATCH --gres=gpu:1
#SBATCH --mem=8g
#SBATCH --time=96:00:00
#SBATCH --partition=ggo

conda activate myenv
srun python3 LTIFM2.py




