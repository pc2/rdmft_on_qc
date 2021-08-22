for s in `seq 1 7`;
do
  for m in JordanWigner BravyiKitaev Parity; 
  do
    for l in disjointqubits qubitwise commute;
    do
      cat combination_of_observables.cfg  | sed "s/_COMMUTATION_/$l/g" | sed "s/_SIZE_/$s/g" | sed "s/_MAPPING_/$m/g" > combination_of_observables_${s}_${m}_${l}.cfg
      python3 ../../rdmft.py combination_of_observables_${s}_${m}_${l}.cfg > combination_of_observables_${s}_${m}_${l}.out
      rm *.png
    done
  done
done
