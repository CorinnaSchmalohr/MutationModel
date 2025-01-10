#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G  
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/WGSanalysisPlots_o-%j
#SBATCH --error=output/WGSanalysisPlots_e-%j
#SBATCH -w beyer-n02


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/05_analysis/01d_modelEvaluationWGS.R' 