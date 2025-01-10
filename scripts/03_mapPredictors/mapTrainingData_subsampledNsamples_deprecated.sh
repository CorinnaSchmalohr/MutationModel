#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/mapPredictorsSubNsamples_o-%j
#SBATCH --error=output/mapPredictorsSubNsamples_e-%j

module load R-4.1.2
Rscript --vanilla --verbose 'scripts/03_mapPredictors/03b_mapTrainingData_subsampledNsamples.R' 
