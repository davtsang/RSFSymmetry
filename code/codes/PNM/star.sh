#!/bin/bash
for files in EoSs/*
do
#head -n -9 $files > test1.txt
#tail -n +2 test1.txt > $files
./TOV_repeat $files 1.40
./star_repeat $files 0 0.1 0.01 5.0 1.0
done
exit 0






