#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --account=nbundsch
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=nbundsch@uni-koeln.de
#SBATCH --output=/cellnet/MutationModel/output/GLM_drop1_Oo-%j
#SBATCH --error=/cellnet/MutationModel/output/GLM_drop1_Eo-%j
#SBATCH --array=1-7

Rscript --vanilla --verbose '/cellnet/MutationModel/scripts/logRegression/04a_logR_cluster_drop1.R' ${SLURM_ARRAY_TASK_ID}
