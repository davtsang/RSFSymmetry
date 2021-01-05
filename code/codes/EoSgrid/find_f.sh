#!/bin/bash
for K in $(seq -200 40 40)
do
for f in $(seq 80 20 300)
do
./find_f $f $K 1.2 2 >"contours_JL/K$K""_M1.4_chirp1.2_dfm_f$f.txt"
./find_f $f $K 1.2 1 >"contours_JL/K$K""_M1.4_chirp1.2_dfp_f$f.txt"
./find_f $f $K 1.2 0 >"contours_JL/K$K""_M1.4_chirp1.2_f$f.txt"
done
done
exit 0


##!/bin/bash
#for K in $(seq 20 10 80)
#do
#for f in $(seq 80 20 300)
#do
#./find_f $f $K 1.2 2 >"contours_JK/L$K""_M1.4_chirp1.2_dfm_f$f.txt"
#./find_f $f $K 1.2 1 >"contours_JK/L$K""_M1.4_chirp1.2_dfp_f$f.txt"
#./find_f $f $K 1.2 0 >"contours_JK/L$K""_M1.4_chirp1.2_f$f.txt"
#done
#done
#exit 0


##!/bin/bash
#for K in $(seq 25 1 35)
#do
#for f in $(seq 80 20 300)
#do
#./find_f $f $K 1.2 2 >"contours_LK/J$K""_M1.4_chirp1.2_dfm_f$f.txt"
#./find_f $f $K 1.2 1 >"contours_LK/J$K""_M1.4_chirp1.2_dfp_f$f.txt"
#./find_f $f $K 1.2 0 >"contours_LK/J$K""_M1.4_chirp1.2_f$f.txt"
#done
#done
#exit 0
