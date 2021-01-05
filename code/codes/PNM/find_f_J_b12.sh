#!/bin/bash
for f in $(seq 100 20 300)
do
#for J in $(seq 25.0 1.4 36.2)
#do
./find_f $f 5.991 0 0 >"contours_J_b12/f$f""_a0_5.991.txt"
./find_f $f 6.107 0 0 >"contours_J_b12/f$f""_a0_6.107.txt"
./find_f $f 6.215 0 0 >"contours_J_b12/f$f""_a0_6.215.txt"
./find_f $f 6.322 0 0 >"contours_J_b12/f$f""_a0_6.322.txt"
./find_f $f 6.438 0 0 >"contours_J_b12/f$f""_a0_6.438.txt"
./find_f $f 6.547 0 0 >"contours_J_b12/f$f""_a0_6.547.txt"
./find_f $f 6.655 0 0 >"contours_J_b12/f$f""_a0_6.655.txt"
./find_f $f 6.761 0 0 >"contours_J_b12/f$f""_a0_6.761.txt"
./find_f $f 6.878 0 0 >"contours_J_b12/f$f""_a0_6.878.txt"

./find_f $f 5.991 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_5.991.txt"
./find_f $f 6.107 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_6.107.txt"
./find_f $f 6.215 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_6.215.txt"
./find_f $f 6.322 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_6.322.txt"
./find_f $f 6.438 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_6.438.txt"
./find_f $f 6.547 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_6.547.txt"
./find_f $f 6.655 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_6.655.txt"
./find_f $f 6.761 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_6.761.txt"
./find_f $f 6.878 1.2 1 >"contours_J_b12/chirp1.2_pdf_f$f""_a0_6.878.txt"

./find_f $f 5.991 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_5.991.txt"
./find_f $f 6.107 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_6.107.txt"
./find_f $f 6.215 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_6.215.txt"
./find_f $f 6.322 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_6.322.txt"
./find_f $f 6.438 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_6.438.txt"
./find_f $f 6.547 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_6.547.txt"
./find_f $f 6.655 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_6.655.txt"
./find_f $f 6.761 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_6.761.txt"
./find_f $f 6.878 1.2 2 >"contours_J_b12/chirp1.2_mdf_f$f""_a0_6.878.txt"
#done
done
exit 0
