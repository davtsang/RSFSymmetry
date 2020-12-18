#!/bin/bash
  cp temp3.txt "all_ctav.txt"
  for J in $(seq 25 1 35)
  do
    for L in $(seq 20 10 80)
    do
      for K in $(seq -200 40 40)
      do
        ./TOV_ctav $J $L $K 1.4
        ./ctav_JLK $J $L $K > temp1.txt
        cat "all_ctav.txt" temp1.txt >temp2.txt
        mv temp2.txt "all_ctav.txt" 
      done
    done
  done
exit 0






