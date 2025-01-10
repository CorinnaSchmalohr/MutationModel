#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=nbundsc1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=nbundsc1@uni-koeln.de
#SBATCH --output=output/glm_sub_o-%j
#SBATCH --error=output/glm_sub_e-%j
#SBATCH --array=1-10


source /data/public/apapada1/Conda/anaconda/bin/activate
conda activate r4-base
Rscript --vanilla --verbose 'scripts/04_trainmodels/04a_glm_subsampled.R' ${SLURM_ARRAY_TASK_ID}
