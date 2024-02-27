#!/bin/bash

#SBATCH --job-name Detection_Rate_MMD
module load R/4.1.3
echo -e '===================================== A ==============================$'
srun -w node15 R -q  --vanilla <  Detection_Rate_MMD.R  &> logs/Detection_Rate_MMD.log &