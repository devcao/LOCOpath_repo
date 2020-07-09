#!/bin/sh

#SBATCH --job-name=pc_adp_ddi_multi_dep_80

#SBATCH --error=pc_adp_hdi_multi_dep_80%j.err

#SBATCH --output=pc_adp_hdi_multi_dep_80%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-10 # run the simulation this many times

module load R/gcc/3.6.0

Rscript  --vanilla ~/hdi_simu/pc_power_simu_multi_80.R 100 80 500 500 multi_pc_Aadp3rd_hdi_simu_dep ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars de-sparse lars_de  ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars Truth lars_Truth  ${SLURM_ARRAY_TASK_ID}
