#!/bin/bash
wd=`pwd`

U="30"
uopt="False"
mueller="False"
dir="U_${U}/Uopt_${uopt}/Mueller_${mueller}/"
rm -rf $dir
mkdir -p $dir

cat template.cfg | sed "s/_U_/$U/g" | sed "s/_uopt_/$uopt/g" | sed "s/_mueller_/$mueller/g" > $dir/input.cfg

cd $dir
python3 -u $wd/../../rdmft.py input.cfg | tee out
cd $wd
