#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=/cellnet/MutationModel/output/histoneENCODE_o-%j
#SBATCH --error=/cellnet/MutationModel/output/histoneENCODE_e-%j
#SBATCH --array=1-10
source /data/public/cschmalo/anaconda3/bin/activate 
conda activate MutModel
# only temporarily removed tissues
tissues=('brain' 'breast' 'colon' 'esophagus' 'kidney' 'liver' 'lung' 'ovary' 'prostate' 'skin')
t=${tissues[$((${SLURM_ARRAY_TASK_ID} - 1))]}
echo $t
bash data/rawdata/histoneENCODE/${t}/UnifyCommand.sh
echo 'done'
