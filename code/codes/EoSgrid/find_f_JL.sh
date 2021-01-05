#!/bin/bash
for K in $(seq -200 40 40)
do
for f in $(seq 70 30 400)
do
./find_f $f $K 1.2 2 >"contours_JL/K$K""_M1.4_chirp1.2_dfm_f$f.txt"
./find_f $f $K 1.2 1 >"contours_JL/K$K""_M1.4_chirp1.2_dfp_f$f.txt"   #note: you must modify find_f.c so that the file is read in as J,L,K,f, and the outputs are J,L,K,f
./find_f $f $K 1.2 0 >"contours_JL/K$K""_M1.4_chirp1.2_f$f.txt"
done
done
exit 0
