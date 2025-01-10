#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/glm_o-%j
#SBATCH --error=output/glm_e-%j
#SBATCH --array=1-10
module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/03a_glm.R' ${SLURM_ARRAY_TASK_ID}
