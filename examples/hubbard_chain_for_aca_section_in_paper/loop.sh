rm result.dat
for l in 0 1 2 3 4 5 6 7 8; do 
  sed "s/_LEVEL_/$l/g" hubbard_chain.cfg > hubbard_chain_$l.cfg
  python3 ../../rdmft.py hubbard_chain_$l.cfg > $l.out;
  Faca=`grep F_aca $l.out | awk 'BEGIN { FS = "=" } ; { print $2 }'`
  Ford=`grep F_reordered $l.out | awk 'BEGIN { FS = "=" } ; { print $2 }'`
  echo "$l $Faca $Ford" >> result.dat
done
