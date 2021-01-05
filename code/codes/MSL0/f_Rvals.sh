#!/bin/bash
for M in $(seq 2.0 -0.05 0.5)
do
./TOV_f_Rvals EoSs/eos_LQDROP_SkXi450_32.0_50.00_-158.73_2.2_Imid_.dat $M
./star_f_Rvals EoSs/eos_LQDROP_SkXi450_32.0_50.00_-158.73_2.2_Imid_.dat 0 0.1 0.01 10.0 1.0
done
exit 0
