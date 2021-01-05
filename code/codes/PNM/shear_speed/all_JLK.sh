#!/bin/bash
  cp temp3.txt "all_ctav.txt"
  for files in ../EoSs/*
  do
    ./TOV_ctav $files 1.4
    ./ctav_JLK $files > temp1.txt
    cat "all_ctav.txt" temp1.txt >temp2.txt
    mv temp2.txt "all_ctav.txt" 
  done
exit 0






