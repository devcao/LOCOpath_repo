#!/bin/sh

#SBATCH --job-name=sis19_screen_50_1000

#SBATCH --error=sis19_screen_50_1000%j.err

#SBATCH --output=sis19_screen_50_1000%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-4 # run the simulation this many times

module load R/gcc/3.6.0


Rscript  --vanilla ~/hdi_simu/rate_sis_screen.R 20 100 200 ${SLURM_ARRAY_TASK_ID}
Rscript  --vanilla ~/hdi_simu/rate_sis_screen.R 50 1000 200 ${SLURM_ARRAY_TASK_ID}
