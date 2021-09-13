for s in `seq 5 5`;
do
  for m in JordanWigner;# BravyiKitaev Parity; 
  do
    for l in commute; #disjointqubits qubitwise commute;
    do
      cat combination_of_observables.cfg  | sed "s/_COMMUTATION_/$l/g" | sed "s/_SIZE_/$s/g" | sed "s/_MAPPING_/$m/g" > combination_of_observables_${s}_${m}_${l}.cfg
      python3 -u ../../rdmft.py combination_of_observables_${s}_${m}_${l}.cfg > combination_of_observables_${s}_${m}_${l}.out
    done
  done
done
