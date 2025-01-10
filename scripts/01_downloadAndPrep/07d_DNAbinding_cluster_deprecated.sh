#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=/cellnet/MutationModel/output/DNAbinding_Unify_o-%j
#SBATCH --error=/cellnet/MutationModel/output/DNAbinding_Unify_e-%j
#SBATCH --array=1-10
source /data/public/cschmalo/anaconda3/bin/activate 
conda activate MutModel
tissues=$(ls data/rawdata/DNAbinding/)
t=$(echo $tissues | cut -d' ' -f ${SLURM_ARRAY_TASK_ID})
echo $t
bash data/rawdata/DNAbinding/${t}/UnifyCommand.sh
echo 'done'
