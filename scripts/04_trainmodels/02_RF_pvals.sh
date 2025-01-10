#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/RF_pv_o-%j
#SBATCH --error=output/RF_pv_e-%j
#SBATCH --array=1-10

Rscript --vanilla --verbose 'scripts/04_trainmodels/02_RF_pvals.R' ${SLURM_ARRAY_TASK_ID}
