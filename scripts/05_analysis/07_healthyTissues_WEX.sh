#!/bin/bash -l
#SBATCH --cpus-per-task=14
#SBATCH --account=cschmal1
#SBATCH --output=output/07_healthyTissues_WEX_e_%j
#SBATCH --error=output/07_healthyTissues_WEX_o_%j


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/05_analysis/07_healthyTissues_WEX.R' 
