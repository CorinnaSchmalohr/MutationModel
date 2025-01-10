#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/exomeMapping_o-%j
#SBATCH --error=output/exomeMapping_e-%j
#SBATCH --array=1-10
module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/03_mapPredictors/04a_wholeExomePositions.R' ${SLURM_ARRAY_TASK_ID}
Rscript --vanilla --verbose 'scripts/03_mapPredictors/04b_wholeExomePositions_postprocessing.R' ${SLURM_ARRAY_TASK_ID}
