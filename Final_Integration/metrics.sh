#!/bin/bash
#SBATCH --job-name=metrics
#SBATCH --account=project_number
#SBATCH --time=40:30:00
#SBATCH --mem=370G
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1
#SBATCH --output=log/error.out
#SBATCH --error=log/error.err

module load python-data
module load cuda
source /projappl/project_2009478/venv/bin/activate
srun python metrics.py
