#!/bin/bash
module load R/4.1.3
echo -e '===================================== A ==============================$'
srun -w node22 R -q  --vanilla <  HIPC_BDRT_maxsize.R  &> logs/HIPC_BDRT_maxsize.log &