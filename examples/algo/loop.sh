
for i in `seq 2 20`;
do
  for f in JordanWigner Parity BravyiKitaev;
  do
    for t in False True;
    do
      echo "$i $f $t"
      cat combination_of_observables.cfg | sed "s/_L_/$i/g" | sed "s/_F_/$f/g" | sed "s/_T_/$t/g" > combination_of_observables_${i}_${f}_${t}.cfg
      python3 ../../rdmft.py combination_of_observables_${i}_${f}_${t}.cfg 2>&1 > combination_of_observables_${i}_${f}_${t}.out
    done
  done
done
