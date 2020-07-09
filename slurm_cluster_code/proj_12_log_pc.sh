#!/bin/sh

#SBATCH --job-name=pc_proj_log_80_dep

#SBATCH --error=pc_proj_log_80_dep%j.err

#SBATCH --output=pc_proj_log_80_dep%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-10 # run the simulation this many times
module load gcc/6.1.0
module load R/gcc/3.6.0

Rscript  --vanilla ~/hdi_simu/log_power_simu_proj_80.R 100 80 500 ${SLURM_ARRAY_TASK_ID}
#Rscript  --vanilla ~/hdi_simu/temp_run.R 100 1000 500 10 
