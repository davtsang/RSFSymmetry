#!/bin/bash
for f in $(seq 100 5 300)
do
cat J_L_K_contours/f$f.txt L_K_J_contours/f$f.txt > planes/$f.txt
done
exit 0
