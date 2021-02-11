#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=/cellnet/MutationModel/output/RFapplication_o-%j
#SBATCH --error=/cellnet/MutationModel/output/RFapplication_e-%j
#SBATCH --array=1-7

Rscript --vanilla --verbose '/cellnet/MutationModel/scripts/RFapplication/RF_application_alltissues_binWise_predictorPvals.R' ${SLURM_ARRAY_TASK_ID}
