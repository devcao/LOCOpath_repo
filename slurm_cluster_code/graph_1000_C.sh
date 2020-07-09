#!/bin/sh

#SBATCH --job-name=C_graph_1K

#SBATCH --error=C_graph_1K%j.err

#SBATCH --output=C_graph_1K%j.out


# SBATCH -p BigMem

#SBATCH -t 1-00:00 # walltime request D-HH:MM

#SBATCH -n 28

#SBATCH -N 1

# SBATCH --array=0-10 # run the simulation this many times

module load R/gcc/3.6.0
module load gcc/6.1.0

Rscript  --vanilla ~/hdi_simu/graph_screen_type_C.R
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars de-sparse lars_de  ${SLURM_ARRAY_TASK_ID}
# Rscript  --vanilla ~/hdi_path/power_simu.R 100 1000 500 500 L2.squared lars Truth lars_Truth  ${SLURM_ARRAY_TASK_ID}
