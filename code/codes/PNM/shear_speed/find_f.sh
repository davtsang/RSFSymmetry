#!/bin/bash
for J in $(seq 25 1 35)
do
for f in $(seq 2000000 500000 10000000)
do
./find_f $f $J >"ctav_contours/$J""_L""_K""_M1.4_ctav$f.txt"
done
done
exit 0

