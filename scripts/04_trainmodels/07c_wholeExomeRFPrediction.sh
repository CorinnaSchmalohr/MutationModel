#!/bin/bash -l
#SBATCH --cpus-per-task=28
#SBATCH --account=cschmal1
#SBATCH --output=output/exomeRF_o-%j
#SBATCH --error=output/exomeRF_e-%j

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/07c_wholeExomeRFPrediction.R'