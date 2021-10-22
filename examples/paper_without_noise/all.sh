#!/bin/bash

if [[ "$1" == "1" ]]; then
  for U in 1; 
  do
    bash run.sh $U True False
    bash run.sh $U True True
    bash run.sh $U False False
  done
  gnuplot plot1.gnu
fi
if [[ "$1" == "2" ]]; then
  #for U in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.5 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
  for U in 18; do #0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.5 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
    bash run.sh $U True True
  done
  bash collect2.sh True True > U.dat
  gnuplot plot2.gnu
fi
