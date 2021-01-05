#!/bin/bash
for f in $(seq 100 5 300)
do
./find_f $f 25.0 > test2.txt
for J in $(seq 26.4 1.4 36.2)
do
cat test2.txt rough_contours/empty.txt >test3.txt
./find_f $f $J > test1.txt
cat test3.txt test1.txt >test2.txt
done
mv test2.txt L_K_J_contours/f$f.txt
done
exit 0
