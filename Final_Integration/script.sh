##Puhti sh scripts submitted to the system, GPU 

#!/bin/bash
#SBATCH --job-name=integratescvigpu
#SBATCH --account=project_2009478
#SBATCH --time=10:30:00
#SBATCH --mem=370G
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1
#SBATCH --output=log/error.out
#SBATCH --error=log/error.err

module load python-data
module load cuda
source /projappl/project_2009478/venv/bin/activate
##replace the name of the script
srun python script.py
