#!/bin/bash
#SBATCH --job-name=Annotation
#SBATCH --account=project_name
#SBATCH --time=00:00:00
#SBATCH --mem=GGG 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task= 00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1
#SBATCH --output=log/test_%j.out
#SBATCH --error=log/test_%j.err

module load python-data
##check the versions of cuda and jax, and make sure you install the compatible versions if needed
module load cuda
module load jax

srun python script.py
