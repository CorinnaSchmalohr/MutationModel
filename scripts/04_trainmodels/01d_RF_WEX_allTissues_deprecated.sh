#!/bin/bash -l
#SBATCH --cpus-per-task=56
#SBATCH --account=cschmalo
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/RF_allTissues_e-%j
#SBATCH --error=output/RF_allTissues_o-%j


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/01d_RF_WEX_allTissues.R' 
