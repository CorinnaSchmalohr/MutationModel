#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=nbundsc1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=nbundsc1@uni-koeln.de
#SBATCH --output=output/RF_o-%j
#SBATCH --error=output/RF_e-%j
#SBATCH --array=1-2


# careful: have to use R4.x
source /data/public/apapada1/Conda/anaconda/bin/activate
conda activate r4-base
Rscript --vanilla --verbose 'scripts/04_trainmodels/09_RF_allTissues.R' ${SLURM_ARRAY_TASK_ID}
