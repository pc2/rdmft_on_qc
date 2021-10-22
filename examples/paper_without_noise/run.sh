#!/bin/bash
wd=`pwd`

U="$1"
uopt="$2"
mueller="$3"
dir="U_${U}/Uopt_${uopt}/Mueller_${mueller}/"
rm -rf $dir
mkdir -p $dir

cat template.cfg | sed "s/_U_/$U/g" | sed "s/_uopt_/$uopt/g" | sed "s/_mueller_/$mueller/g" > $dir/input.cfg

cd $dir
python3 -u $wd/../../rdmft.py input.cfg | tee out

Faca=`grep F_aca out | awk 'BEGIN { FS = " " } ; { print $5 }' `
Freordered=`grep F_reordered out | awk 'BEGIN { FS = " " } ; { print $5 }' `
echo "#Faca=$Faca Fcluster=$Freordered" > conv
echo "#L W sum(c^2)" >> conv
grep -B 1 point out  | grep -v point | grep -v "\-\-" | awk 'BEGIN { FS = " " } ; { print $4,$6,$8 }' >> conv
cd $wd
