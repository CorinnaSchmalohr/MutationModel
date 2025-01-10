#!/bin/bash -l
#SBATCH --cpus-per-task=48
#SBATCH --account=cschmalo
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/RFWGS_o-%j
#SBATCH --error=output/RFGWS_e-%j
#SBATCH --array=8


module unload R-3.5.1
module load R-4.1.2
# conda activate Rbase
Rscript --vanilla --verbose 'scripts/04_trainmodels/01b_RF_WGSdata_skin.R' ${SLURM_ARRAY_TASK_ID}
