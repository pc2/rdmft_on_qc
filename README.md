# RDMF_on_QC

[![DOI](https://zenodo.org/badge/434230783.svg)](https://zenodo.org/badge/latestdoi/434230783)


*Reduced density matrix functional evaluation on a quantum computer.*

## Installation
### Linux

Make sure you have Python (version >= 3.8) installed. Then simply run `sh env.sh` to create a new virtual environment and install all necessary depencencies.

### MacOS

First follow the Linux instructions above. Note, however, that the C++ code in dmrgpy will fail to compile. In the following, we'll fix this.

You will need to compile `deps/dmrgpy/src/mpscpp2/ITensor` manually. Given that a GNU g++ compiler is available (can be installed via `brew install gcc`), copy the file `options.save` to `options.mk` and make sure to set `CCCOM=g++-11 -m64 -std=c++14 -fPIC` (comment out the previous linux line). Note that using clang will work for ITensor but not for the next step.

After we've compiled ITensor, we go up one directory (to `deps/dmrgpy/src/mpscpp2`) and build this as well: `make -j 4`. After completion execute `mv mpscpp mpscpp.x`. Finally, make sure to set `export DMRGROOT=/path/to/repo/deps/dmrgpy/src` either in your `.bashrc` or (to keep things contained) in a local `.envrc` in combination with using [direnv](https://direnv.net/).

That's it, you should be good to go.
## Usage

The main python script is `rdmft.py` which accepts a configuration file (`.cfg`) like `hubbard_chain.cfg`. Simply run

```
python3 rdmft.py hubbard_chain.cfg
```

## What does the script do?
1. sets up a system like a Hubbard chain
2. computes the numerically exact ground state with MPS and the corresponding one-particle reduced density matrix (1RDM)
3. applies the local approximation to the density-matrix functional (Eq. 5.34 in https://ediss.uni-goettingen.de/bitstream/handle/11858/00-1735-0000-002E-E5C2-7/out.pdf)
4. applies the adaptive cluster approximation (ACA) to the ilocal-th local density-matrix functional (chapter 8 in https://ediss.uni-goettingen.de/bitstream/handle/11858/00-1735-0000-002E-E5C2-7/out.pdf)
5. computes the numerically exact value of the ilocal-th local density-matrix functional in the ACA with a full state-vector parametrization (aka. ED, aka. FCI) of the wave function
6. sets up a parametrized hardware-efficient trial state as a variational state
7. sets up all quantum programs including measurements for the required observables
8. computes the ilocal-th local density-matrix functional in the ACA by constrained minimization with the augmented Lagrangian over the parameters of the hardware-efficient trial state (get multiple coffees)
  * Either with
    * simulation of the quantum programs on a classical computer without sampling, without noise and with explicit constraints (tsim=True, tsampling=False, tnoise=False, tsimcons=True)
    * simulation of the quantum programs on a classical computer with sampling, without noise and with implicit constraints (tsim=True, tsampling=True, tnoise=False, tsimcons=False)
      * currently without convergence criterion fou outer iterations of augmented Lagrangian
    * simulation of the quantum programs on a classical computer with sampling, with noise and with implicit constraints (tsim=True, tsampling=True, tnoise=True, tsimcons=False)
      * currently without convergence criterion fou outer iterations of augmented Lagrangian
    * simulation of the quantum programs on a quantum computer (tsim=False, tsampling=True, tnoise=True, tsimcons=False)
      * currently without convergence criterion fou outer iterations of augmented Lagrangian

