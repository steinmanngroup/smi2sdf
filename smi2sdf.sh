#!/usr/bin/env bash
#SBATCH --job-name rdkitconf
#SBATCH --partition teach
#SBATCH --time 24:00:00
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --array 1-260

python smi2sdf.py ZINC_250k.smi ${SLURM_ARRAY_TASK_ID} 1000
