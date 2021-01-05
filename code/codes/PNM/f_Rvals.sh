#!/bin/bash
for M in $(seq 2.0 -0.05 0.5)
do
./TOV_f_Rvals 31.0 40.01 13 $M
./star_f_Rvals 31.0 40.01 13 0 0.2 0.01 20.0 1.0
done
exit(0)
