#!/bin/bash
module load R/4.1.3
echo -e '===================================== A ==============================$'
srun -w node15 R -q  --vanilla <  HIPC_MMD_maxsize_9.R  &> logs/HIPC_MMD_maxsize_9.log &