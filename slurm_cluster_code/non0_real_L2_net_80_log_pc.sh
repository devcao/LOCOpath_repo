#!/bin/sh

#SBATCH --job-name=pcL2_non0_8log_80_dep

#SBATCH --error=pcL2_non0_8log_80_dep%j.err

#SBATCH --output=pcL2_non0_8log_80_dep%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-10 # run the simulation this many times

module load gcc/6.1.0
module load R/gcc/3.6.0

Rscript  --vanilla ~/hdi_simu/non0_pc_net_power_simu_80_log_L2.R 100 80 500 500 non0_log_pc_net_hdi_simu_80_dep ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars de-sparse lars_de  ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars Truth lars_Truth  ${SLURM_ARRAY_TASK_ID}
