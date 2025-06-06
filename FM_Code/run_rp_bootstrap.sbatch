#!/bin/bash --login
#SBATCH --job-name=bootstrap_validation
#SBATCH --time=96:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=48

conda activate venv
srun python3 rp_bootstrap_mp.py




