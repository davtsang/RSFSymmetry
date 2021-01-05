#!/bin/bash
#                                                                                              HOW TO MODIFY TO GET ONE OF THE OTHER PARAMETERS AS THE CONSTANT ONE:
  for J in $(seq 25 1 35)                                                                    #non-constant parameter loop
  do
    echo $J                                                                                  #change to the non-constant parameter
    cp temp3.txt "$J""_L""_ctav_M16.txt"                                                      #put $ in front of the non-constant parameter, and not the constant one
    for L in $(seq 20 10 90)                                                                 #constant parameter loop
    do
      K=$(($((371*$L))-$((1113*$J))+1193))
#      echo $K
      ./TOV_ctav $J $L $K 1.6
      ./ctav_JLK $J $L $K > temp1.txt
      cat "$J""_L""_ctav_M16.txt" temp1.txt >temp2.txt                                        #put $ in front of the non-constant parameters, and not the constant one
      mv temp2.txt "$J""_L""_ctav_M16.txt"                                                    #put $ in front of the non-constant parameters, and not the constant one
    done
  done
exit 0






