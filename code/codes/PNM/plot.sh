#!/bin/bash
for f in $(seq 100 5 300)
do
cat "rough_contours/f$f""_K-300.txt" "rough_contours/f$f""_K-295.txt" >test1.txt
for K in $(seq -290 5 0)
do
cat test1.txt "rough_contours/f$f""_K$K.txt" >test2.txt
mv test2.txt test1.txt
done
mv test1.txt J_L_K_contours/f$f.txt
done
exit 0
