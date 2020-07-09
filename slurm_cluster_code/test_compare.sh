#!/bin/sh

#SBATCH --job-name=run_compare

#SBATCH --error=run_comp%j.err

#SBATCH --output=run_comp%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0 # run the simulation this many times

module load R/gcc/3.6.0

Rscript  --vanilla ~/hdi_simu/run_compare.R 
