#!/bin/bash
#                                                                                                HOW TO MODIFY TO GET ONE OF THE OTHER PARAMETERS AS THE CONSTANT ONE:
  for J in $(seq 25 1 35)                                                                                         #non-constant parameter loop 1
  do
    for L in $(seq 20 10 80)                                                                                      #non-constant parameter loop 2
    do
      echo $J $L                                                                                                  #change to the two non-constant parameters
      cp temp3.txt "$J""_$L""_K""_ctav.txt"                                                                       #put $ in front of the non-constant parameters, and not the constant one
      for K in $(seq -200 40 40)                                                                                  #constant parameter loop
      do
        ./TOV_ctav $J $L $K 1.4                                                                                   #change to the non-constant parameters
        ./ctav_JLK $J $L $K > temp1.txt
        cat "$J""_$L""_K""_ctav.txt" temp1.txt >temp2.txt                                                         #put $ in front of the non-constant parameters, and not the constant one
        mv temp2.txt "$J""_$L""_K""_ctav.txt"                                                                     #put $ in front of the non-constant parameters, and not the constant one
      done
    done
  done
exit 0






