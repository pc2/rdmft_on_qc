# RDMF_on_QC

# Dependencies

* Python >= 3.8
* https://github.com/joselado/dmrgpy: `python install.py`
* https://github.com/quantumlib/OpenFermion: `python -m pip install --user openfermion`
* https://github.com/Qiskit/qiskit-nature: `python -m pip install --user qiskit-nature`

# Setup of Depencencies
* simply run `bash env.sh`

# What This Script Does
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

# What is Currently  Missing:
* Pulse mode https://git.uni-paderborn.de/pc2/quantum-computing/nhr-qc/rdmf_on_qc/-/issues/10
* hardware-efficient placement of cnots in measurements https://git.uni-paderborn.de/pc2/quantum-computing/nhr-qc/rdmf_on_qc/-/issues/13: partially done
  * for measurement circuits that measure only a single qubit
  * implicitly for construction of measurement circuits of commuting measurements
* testing

# How to Use
```
python3 rdmft.py hubbard_chain.cfg
```
