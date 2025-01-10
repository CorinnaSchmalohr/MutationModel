#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=nbundsc1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=nbundsc1@uni-koeln.de
#SBATCH --output=output/chr1_allMuts_o-%j
#SBATCH --error=output/chr1_allMuts_e-%j
#SBATCH --array=1-10

source /data/public/apapada1/Conda/anaconda/bin/activate
conda activate r4-base
conda install -c r r-xlsx
Rscript --vanilla --verbose 'scripts/03_mapPredictors/chr1arm_allMuts_cluster.R' ${SLURM_ARRAY_TASK_ID}