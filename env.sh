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

cd $wd


echo "source $envdir/bin/activate"


