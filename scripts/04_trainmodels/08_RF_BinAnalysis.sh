#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=nbundsc1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=nbundsc1@uni-koeln.de
#SBATCH --output=/cellnet/MutationModel/output/RF_Bin_o-%j
#SBATCH --error=/cellnet/MutationModel/output/RF_Bin_e-%j
#SBATCH --array=1-2

source /data/public/apapada1/Conda/anaconda/bin/activate
conda activate r4-base
Rscript --vanilla --verbose '/cellnet/MutationModel/scripts/04_trainmodels/08_RF_BinAnalysis_bash.R' ${SLURM_ARRAY_TASK_ID}
