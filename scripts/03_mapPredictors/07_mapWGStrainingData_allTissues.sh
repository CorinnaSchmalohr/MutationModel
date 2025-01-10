#!/bin/bash -l
#SBATCH --account=cschmal1
#SBATCH --mem 250G
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/mapAllTissueWGSmuts_o-%j
#SBATCH --error=output/mapAllTissueWGSmuts_e-%j
#SBATCH --partition=all

module unload R-3.5.1 
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/03_mapPredictors/07_mapWGStrainingData_allTissues.R' 
