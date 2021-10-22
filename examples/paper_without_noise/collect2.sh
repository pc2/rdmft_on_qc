uopt="$1"
mueller="$2"

for U in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.5 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
  dir="U_${U}/Uopt_${uopt}/Mueller_${mueller}/"
  Faca=`grep "Faca" $dir/conv | awk 'BEGIN { FS = "=" } ; { print $2 }' | awk 'BEGIN { FS = " " } ; { print $1 }'`
  F=`tail -n 1 $dir/conv | awk 'BEGIN { FS = " " } ; { print $2 }'`
  c=`tail -n 1 $dir/conv | awk 'BEGIN { FS = " " } ; { print $3 }'`
  echo "$U $Faca $F $c"
done
