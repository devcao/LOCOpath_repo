#!/bin/sh

#SBATCH --job-name=poi_Tr_10net_hdi_1000_dep

#SBATCH --error=poi_Tr_10net_hdi_1000_dep%j.err

#SBATCH --output=poi_Tr_10net_hdi_1000_dep%j.out


#SBATCH -p BigMem

#SBATCH -t 15-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

#SBATCH --array=0 # run the simulation this many times

module load gcc/6.1.0
module load R/gcc/3.6.0

Rscript  --vanilla ~/hdi_simu/pc_net_power_simu_1000_10_poi.R 100 1000 500 500 poi_Truth_new_pc_net_hdi_simu_1000_dep ${SLURM_ARRAY_TASK_ID}

# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars de-sparse lars_de  ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars Truth lars_Truth  ${SLURM_ARRAY_TASK_ID}
