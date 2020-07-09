#!/bin/sh

#SBATCH --job-name=var_screen_70_1000

#SBATCH --error=var_screen_70_1000%j.err

#SBATCH --output=var_screen_70_1000%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-3 # run the simulation this many times

module load R/gcc/3.6.0


Rscript  --vanilla ~/hdi_simu/rate_screen.R 70 1000 200 ${SLURM_ARRAY_TASK_ID}
