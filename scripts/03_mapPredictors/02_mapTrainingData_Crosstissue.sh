#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/mapPredictorsCrossTissue_o-%j
#SBATCH --error=output/mapPredictorsCrossTissue_e-%j
#SBATCH --array=1-10

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/03_mapPredictors/02_mapTrainingData_Crosstissue.R'  ${SLURM_ARRAY_TASK_ID}
