#!/bin/bash
module load R/4.1.3
echo -e '===================================== A ==============================$'
srun -w node22 R -q  --vanilla <  DROP_ABC.R  &> logs/DROP_ABC.log &