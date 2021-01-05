#!/bin/bash
for M in $(seq 2.0 -0.05 0.8)
do
./TOV_fR EoSs/eos_LQDROP_SkXi450_30.0_50.00_-120.00_2.2_Imid_.dat $M
./star_fR EoSs/eos_LQDROP_SkXi450_30.0_50.00_-120.00_2.2_Imid_.dat 0 0.1 0.01 10.0 1.0
done
exit 0
