#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/chr1_o-%j
#SBATCH --error=output/chr1_e-%j
#SBATCH --array=1-10


source /data/public/apapada1/Conda/anaconda/bin/activate
conda activate r4-base
Rscript --vanilla --verbose 'scripts/05_analysis/chr1_predictions.R' ${SLURM_ARRAY_TASK_ID}
