#!/bin/bash -l
#SBATCH --cpus-per-task=14
#SBATCH --account=cschmalo
#SBATCH --output=output/09a_generalModelEvaluation_WGS_o_%j
#SBATCH --error=output/09a_generalModelEvaluation_WGS_e_%j

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose "scripts/05_analysis/09b_generalModelEvaluation_WGS.R"
