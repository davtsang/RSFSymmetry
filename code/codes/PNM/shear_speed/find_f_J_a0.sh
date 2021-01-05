#!/bin/bash
for f in $(seq 4000000 250000 7000000)
do
#for J in $(seq 25.0 1.4 36.2)
#do
./find_f $f 5.755 >"ctav_contours/b12_5.755_ctav$f.txt"
./find_f $f 7.786 >"ctav_contours/b12_7.786_ctav$f.txt"
./find_f $f 9.815 >"ctav_contours/b12_9.815_ctav$f.txt"
./find_f $f 11.80 >"ctav_contours/b12_11.80_ctav$f.txt"
./find_f $f 13.79 >"ctav_contours/b12_13.79_ctav$f.txt"
./find_f $f 15.80 >"ctav_contours/b12_15.80_ctav$f.txt"
./find_f $f 17.80 >"ctav_contours/b12_17.80_ctav$f.txt"
./find_f $f 19.76 >"ctav_contours/b12_19.76_ctav$f.txt"
./find_f $f 21.81 >"ctav_contours/b12_21.81_ctav$f.txt"
#done
done
exit 0
