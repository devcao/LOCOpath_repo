#!/bin/sh

#SBATCH --job-name=run_shrt_ribo

#SBATCH --error=run_shrt_ribo%j.err

#SBATCH --output=run_shrt_ribo%j.out


#SBATCH -p BigMem

#SBATCH -t 2-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0 # run the simulation this many times

module load R/gcc/3.6.0

Rscript  --vanilla ~/hdi_simu/run_ribo.R 
