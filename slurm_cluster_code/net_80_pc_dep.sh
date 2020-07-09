#!/bin/sh

#SBATCH --job-name=pc_net_hdi_80_dep

#SBATCH --error=pc_net_hdi_80_dep%j.err

#SBATCH --output=pc_net_hdi_80_dep%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-10 # run the simulation this many times

module load R/gcc/3.6.0
module load gcc/6.1.0

Rscript  --vanilla ~/hdi_simu/pc_net_power_simu_80_dep.R 100 80 500 500 ext_pc_net_Aadp_hdi_simu_80_dep ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars de-sparse lars_de  ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars Truth lars_Truth  ${SLURM_ARRAY_TASK_ID}
