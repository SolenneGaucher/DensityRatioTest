#!/bin/bash
module load R/4.1.3
echo -e '===================================== A ==============================$'
srun -w node15 R -q  --vanilla <  Robustness_test.R  &> logs/Robustness_test.log &