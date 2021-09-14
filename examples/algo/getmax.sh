
for i in `seq 2 10`;
do
  m=`grep transpiled_complexity combination_of_observables_${i}_$1.out | awk 'BEGIN { FS = " ";m=0; } ; {if($7>m){m=$7}} END { print m }'`
  echo "$i $m"

done
