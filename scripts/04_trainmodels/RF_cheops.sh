#!/bin/bash -l
#SBATCH --cpus-per-task=12
#SBATCH --mem=46gb
#SBATCH --time=48:00:00
#SBATCH --account=cschmal1
#SBATCH --mail-type=ALL
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/RF_o-%j
#SBATCH --error=output/RF_e-%j
#SBATCH --array=1-10

export R_LIBS_USER=$HOME/R/4.0.2
module load R/4.0.2_gnu_mkl
Rscript --vanilla --verbose 'scripts/01_RF.R' ${SLURM_ARRAY_TASK_ID}
