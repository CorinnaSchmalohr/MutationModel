#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/RF_o-%j
#SBATCH --error=output/RF_e-%j
#SBATCH --array=1-10


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/01a_RF.R' ${SLURM_ARRAY_TASK_ID}
