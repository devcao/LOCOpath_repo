#!/bin/sh

#SBATCH --job-name=eq12_hdi_80_dep

#SBATCH --error=eq12_hdi_80_dep%j.err

#SBATCH --output=eq12_hdi_80_dep%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-10 # run the simulation this many times

module load R/gcc/3.6.0
module load gcc/6.1.0

Rscript  --vanilla ~/hdi_simu/pc_eq_power_simu_12.R 100 12 500 500 pc_right_eq_hdi_simu_12_dep ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars de-sparse lars_de  ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars Truth lars_Truth  ${SLURM_ARRAY_TASK_ID}
