#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=jkrueg17
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=jkrueg17@smail.uni-koeln.de
#SBATCH --output=output/prepareWGSmuts_o-%j
#SBATCH --error=output/prepareWGSmuts_e-%j
#SBATCH --array=1-6

module unload R-3.5.1
module load R-4.1.2
# source /data/public/cschmalo/anaconda3/bin/activate 
# conda activate test
Rscript --vanilla --verbose 'scripts/02_prepareMutations/06b_prepare_scWGSforModelevaluation_new.R'  ${SLURM_ARRAY_TASK_ID}
