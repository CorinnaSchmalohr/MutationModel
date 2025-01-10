#!/bin/bash -l
#SBATCH --account=cschmal1
#SBATCH --cpus-per-task=56
#SBATCH --mem 250G
#SBATCH --output=output/01_mapTrainingData_allTissues_o_%j
#SBATCH --error=output/01_mapTrainingData_allTissues_e_%j
#SBATCH --partition=all

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/03_mapPredictors/01_mapTrainingData_allTissues.R'
