#!/bin/bash
for K in $(seq 20 10 80)
do
for f in $(seq 70 30 400)
do
./find_f $f $K 1.2 2 >"contours_JK/L$K""_M1.4_chirp1.2_dfm_f$f.txt"
./find_f $f $K 1.2 1 >"contours_JK/L$K""_M1.4_chirp1.2_dfp_f$f.txt"   #note: you must modify find_f.c so that the file is read in as J,K,L,f, and the outputs are J,K,L,f
./find_f $f $K 1.2 0 >"contours_JK/L$K""_M1.4_chirp1.2_f$f.txt"
done
done
exit 0
