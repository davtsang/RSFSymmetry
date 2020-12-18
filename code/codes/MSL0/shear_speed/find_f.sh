#!/bin/bash
  for f in $(seq 40000000 2500000 100000000)
  do
    ./find_f $f 0 >"ct_contours_nb37/ct$f.txt"
  done
exit 0

