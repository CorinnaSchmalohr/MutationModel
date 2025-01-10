#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=nbundsc1@uni-koeln.de
#SBATCH --output=output/RFsubsampledNsamples_o-%j
#SBATCH --error=output/RFsubsampledNsamples_e-%j
#SBATCH --array=1-10


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/04_RF_subsampledNsamples.R' ${SLURM_ARRAY_TASK_ID}

