export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

wd=`pwd`
for i in  `seq 2 24`;
do
  for f in JordanWigner Parity BravyiKitaev;
  do
    for t in False;
    do
      echo "$i $f $t"
      cd $wd
      dir="combination_of_observables_${i}_${f}_${t}"
      mkdir $dir
      cat combination_of_observables.cfg | sed "s/_L_/$i/g" | sed "s/_F_/$f/g" | sed "s/_T_/$t/g" > $dir/inputf.cfg
      cd $dir
      python3 -u $wd/../../rdmft.py inputf.cfg $i $f $t >out &
      l=`ps aux | grep inputf | wc -l`
      while  [ $l -gt 40 ] ;do
        sleep 1
        echo $l
        l=`ps aux | grep inputf | wc -l`
      done
    done
  done
done
