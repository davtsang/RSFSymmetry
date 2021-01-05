#!/bin/bash
for f in $(seq 4000000 250000 7000000)
do
./find_f $f 5.991 >"ctav_contours/a0_5.991_ctav$f.txt"
./find_f $f 6.107 >"ctav_contours/a0_6.107_ctav$f.txt"
./find_f $f 6.215 >"ctav_contours/a0_6.215_ctav$f.txt"
./find_f $f 6.322 >"ctav_contours/a0_6.322_ctav$f.txt"
./find_f $f 6.438 >"ctav_contours/a0_6.438_ctav$f.txt"
./find_f $f 6.547 >"ctav_contours/a0_6.547_ctav$f.txt"
./find_f $f 6.655 >"ctav_contours/a0_6.655_ctav$f.txt"
./find_f $f 6.761 >"ctav_contours/a0_6.761_ctav$f.txt"
./find_f $f 6.878 >"ctav_contours/a0_6.878_ctav$f.txt"
done
exit 0
