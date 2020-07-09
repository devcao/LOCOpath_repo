#!/bin/sh

#SBATCH --job-name=poi_M3sig_hdi_multi_dep

#SBATCH --error=poi_M3sig_hdi_multi_dep%j.err

#SBATCH --output=poi_M3sig_hdi_multi_dep%j.out


# SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0-10 # run the simulation this many times

module load R/gcc/3.6.0
module load gcc/3.6.0

Rscript  --vanilla ~/hdi_simu/pc_power_simu_multi_poi.R 100 1000 500 500 poi_poi_3signal_multi_pc_hdi_simu_dep ${SLURM_ARRAY_TASK_ID}

# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars de-sparse lars_de  ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars Truth lars_Truth  ${SLURM_ARRAY_TASK_ID}
