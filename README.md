# RDMF_on_QC

# Dependencies

* https://github.com/joselado/dmrgpy: python install.py 
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
    * simulation of the quantum programs on a classical computer without noise (tsim=True, tnoise=False)
    * simulation of the quantum programs on a classical computer with noise (tsim=True, tnoise=True)
    * execution of the quantum programs on a quantum computer (tsim=False))

# What is Currently  Missing:
* possible reuse of Pauli-measurements for multiple observables https://git.uni-paderborn.de/pc2/quantum-computing/nhr-qc/rdmf_on_qc/-/issues/12
* reduction of number of quantum programs by combining multiple observables into one program https://git.uni-paderborn.de/pc2/quantum-computing/nhr-qc/rdmf_on_qc/-/issues/11
* Pulse mode https://git.uni-paderborn.de/pc2/quantum-computing/nhr-qc/rdmf_on_qc/-/issues/10
* hardware-efficient placement of cnots in measurements https://git.uni-paderborn.de/pc2/quantum-computing/nhr-qc/rdmf_on_qc/-/issues/13
* testing

# How to Use
```
python3 rdmft.py hubbard_chain.cfg
```
