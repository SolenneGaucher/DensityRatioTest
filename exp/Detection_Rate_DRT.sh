#!/bin/bash
module load R/4.1.3
echo -e '===================================== A ==============================$'
srun -w node13 R -q  --vanilla <  Detection_Rate_DRT.R  &> logs/Detection_Rate_DRT.log &