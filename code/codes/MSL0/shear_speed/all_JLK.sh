#!/bin/bash
  cp temp3.txt "all_ct_nb37.txt"
  for J in $(seq 25 1 35)
  do
    for L in $(seq 20 10 90)
    do
      K=$(($((371*$L))-$((1113*$J))+1193))
      ./TOV_ctav $J $L $K 1.4                                                                       #change to the non-constant parameters      #in TOV_mucc.c : change letter
      ./ctav_JLK $J $L $K > temp1.txt  
      cat "all_ct_nb37.txt" temp1.txt >temp2.txt
      mv temp2.txt "all_ct_nb37.txt" 
    done
  done
exit 0
