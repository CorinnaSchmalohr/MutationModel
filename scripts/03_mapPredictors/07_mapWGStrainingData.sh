#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/mapWGSPredictors_o-%j
#SBATCH --error=output/mapWGSPredictors_e-%j
#SBATCH --array=8
# 1-
module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/03_mapPredictors/07_mapWGStrainingData.R' ${SLURM_ARRAY_TASK_ID}
