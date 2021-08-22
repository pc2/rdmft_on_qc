PYCMD="python3"
PIPCMD="python3 -m pip"
wd=`pwd`
envdir="$wd/rdmft-qc-env"
depdir="$wd/deps"

#deactivate if we are in a python env
deactivate

#cleanup existing env and create new
if [[ -d "$envdir" ]]; then
  echo "Remove existing env at $envdir and start new? [y/N]"
  read a
  if [[ "$a" == "y" ]];
  then
    rm -rf $envdir
    $PYCMD -m venv $envdir
  fi
else 
  $PYCMD -m venv $envdir
fi

source $envdir/bin/activate


if [[ -d "$depdir" ]]; then
  echo "Remove existing dependencies at $depdir and start new? [y/N]"
  read a
  if [[ "$a" == "y" ]];
  then
    rm -rf $depdir
    mkdir $depdir
  fi
else
  mkdir $depdir
fi

#install packages
if [[ ! -e "requirements.txt" ]]; then
  echo "installing fresh"
  $PIPCMD install wheel
  $PIPCMD install numpy
  $PIPCMD install scipy
  $PIPCMD install python-dotenv
  $PIPCMD install matplotlib
  $PIPCMD install qiskit
  $PIPCMD install qiskit-nature
  $PIPCMD install pylatexenc
  $PIPCMD install sympy
  $PIPCMD install configparser
  $PIPCMD install networkx
  $PIPCMD install numba
else
  echo "installing according to requirements.txt"
  $PIPCMD install -r requirements.txt
fi

#install dmrgpy
if [[ ! -d "$depdir/dmrgpy" ]];then
  cd $depdir
  git clone https://github.com/joselado/dmrgpy
  cd dmrgpy
  git checkout 8895a2ae4f613cdddeb8501fee5eae999e5cb84c
  $PYCMD install.py
fi

#install OpenFermion
if [[ ! -d "$depdir/OpenFermion" ]];then
  cd $depdir
  git clone https://github.com/quantumlib/OpenFermion
  cd OpenFermion
  git checkout 84d169b529a27a0cfd1cc41073d2665f1fcf0312
  $PYCMD -m pip install openfermion
fi

#install vqe-term-grouping
if [[ ! -d "$depdir/vqe-term-grouping" ]];then
  cd $depdir
  #git clone https://github.com/teaguetomesh/vqe-term-groupinig
  git clone https://git.uni-paderborn.de/pc2/quantum-computing/nhr-qc/vqe-term-grouping.git
  cd vqe-term-grouping
  #git checkout 48abb10122a9cd861cc9f3ef3683c4c03b100a77
  git checkout ead475bb9b9a130407bd89c7624c89b9d212e1fd
fi


cd $wd


echo "source $envdir/bin/activate"


