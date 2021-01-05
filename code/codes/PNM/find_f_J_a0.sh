#!/bin/bash
for f in $(seq 100 20 300)
do
#for J in $(seq 25.0 1.4 36.2)
#do
./find_f $f 5.755 0 0 >"contours_J_a0/f$f""_b12_5.755.txt"
./find_f $f 7.786 0 0 >"contours_J_a0/f$f""_b12_7.786.txt"
./find_f $f 9.815 0 0 >"contours_J_a0/f$f""_b12_9.815.txt"
./find_f $f 11.80 0 0 >"contours_J_a0/f$f""_b12_11.80.txt"
./find_f $f 13.79 0 0 >"contours_J_a0/f$f""_b12_13.79.txt"
./find_f $f 15.80 0 0 >"contours_J_a0/f$f""_b12_15.80.txt"
./find_f $f 17.80 0 0 >"contours_J_a0/f$f""_b12_17.80.txt"
./find_f $f 19.76 0 0 >"contours_J_a0/f$f""_b12_19.76.txt"
./find_f $f 21.81 0 0 >"contours_J_a0/f$f""_b12_21.81.txt"

./find_f $f 5.755 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_5.755.txt"
./find_f $f 7.786 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_7.786.txt"
./find_f $f 9.815 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_9.815.txt"
./find_f $f 11.80 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_11.80.txt"
./find_f $f 13.79 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_13.79.txt"
./find_f $f 15.80 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_15.80.txt"
./find_f $f 17.80 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_17.80.txt"
./find_f $f 19.76 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_19.76.txt"
./find_f $f 21.81 1.2 1 >"contours_J_a0/chirp1.2_pdf_f$f""_b12_21.81.txt"

./find_f $f 5.755 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_5.755.txt"
./find_f $f 7.786 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_7.786.txt"
./find_f $f 9.815 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_9.815.txt"
./find_f $f 11.80 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_11.80.txt"
./find_f $f 13.79 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_13.79.txt"
./find_f $f 15.80 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_15.80.txt"
./find_f $f 17.80 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_17.80.txt"
./find_f $f 19.76 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_19.76.txt"
./find_f $f 21.81 1.2 2 >"contours_J_a0/chirp1.2_mdf_f$f""_b12_21.81.txt"
#done
done
exit 0
