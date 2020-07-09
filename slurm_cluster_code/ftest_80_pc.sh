#!/bin/sh

#SBATCH --job-name=pc_ftest_80_dep

#SBATCH --error=pc_Fhdi3rd_ftest_dep_80%j.err

#SBATCH --output=pc_Fhdi3rd_ftest_dep_80%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-10 # run the simulation this many times

module load R/gcc/3.6.0

Rscript  --vanilla ~/hdi_simu/power_simu_fftest_80.R 100 80 500 ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars de-sparse lars_de  ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars Truth lars_Truth  ${SLURM_ARRAY_TASK_ID}
