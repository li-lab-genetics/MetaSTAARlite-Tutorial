#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-00:00
#SBATCH --mem=20000
#SBATCH --mail-type=NONE
#SBATCH --array=1-2

module purge
module load r/4.3.1
export R_LIBS_USER=/proj/xihaoli_lab/users/xihaoli/R-4.3.1
R --vanilla --args ${SLURM_ARRAY_TASK_ID} < $1 > "${1}.${SLURM_ARRAY_TASK_ID}.out" 2> "${1}.${SLURM_ARRAY_TASK_ID}.err"

