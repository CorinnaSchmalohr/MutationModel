#!/bin/bash -l
#SBATCH --cpus-per-task=14
#SBATCH --account=cschmalo
#SBATCH --output=output/09a_generalModelEvaluation_WEX_o_%j
#SBATCH --error=output/09a_generalModelEvaluation_WEX_e_%j

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose "scripts/05_analysis/09a_generalModelEvaluation_WEX.R"
