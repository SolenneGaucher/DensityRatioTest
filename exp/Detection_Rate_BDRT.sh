#!/bin/bash

module load R/4.1.3
echo -e '===================================== A ==============================$'
srun -w node22 R -q  --vanilla <  Detection_Rate_BDRT.R  &> logs/Detection_Rate_BDRT.log &