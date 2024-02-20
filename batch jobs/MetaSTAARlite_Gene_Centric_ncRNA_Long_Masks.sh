#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-00:00
#SBATCH --mem=10000
#SBATCH --mail-type=NONE

module purge
module load r/4.3.1
export R_LIBS_USER=/proj/xihaoli_lab/users/xihaoli/R-4.3.1
R --vanilla --args ${SLURM_JOB_ID} < $1 > "${1}.${SLURM_JOB_ID}.out" 2> "${1}.${SLURM_JOB_ID}.err"

