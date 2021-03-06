import math
from dmrgpy import fermionchain
import numpy as np

def hubbard_chain(n,t,U,mu,Lint,mode="ED"):
    """
    Constructs and solves a Hubbard chain with open
    boundary conditions using ED or DMRG.

    Inputs:
    `n`: # of sites, i.e. length of the chain
    `t`: nn-hopping amplitude
    `U`: interaction amplitude
    `mu`: chemical potential
    `Lint`: first `Lint` sites will host a Hubbard interaction
    `mode` (optional): either "ED" or "DMRG"

    Output:
    `[E,D,W]` where
    `E`: groundstate energy
    `D`: one-particle reduced density matrix
    `W`: double occupancy
    """
    fc = fermionchain.Spinful_Fermionic_Chain(n)

    if n==2 and Lint==2:
        #dimer comparison results
        E=-2*t*(math.sqrt(1+(U/(4*t))*(U/(4*t)))-U/(4*t))
        theta=math.atan(math.sqrt(1+(U/(4*t))*(U/(4*t)))+U/(4*t))-math.pi/4
        rho12=0.5*math.cos(2*theta)
        print("E_dimer=",E)
        print("rho12_dimer=",rho12)

    # first neighbor hopping
    h = 0
    if n>1:
        for i in range(n-1):
            h = h + t*fc.Cdagup[i]*fc.Cup[i+1]
            h = h + t*fc.Cdagdn[i]*fc.Cdn[i+1]
        h = h + h.get_dagger() # Make Hermitian
    # Hubbard term
    for i in range(Lint):
      h = h + U*fc.Nup[i]*fc.Ndn[i]

    # chemical potential
    for i in range(n):
      h =h + mu*(fc.Nup[i]+fc.Ndn[i])

    fc.set_hamiltonian(h) # initialize the Hamiltonian
    fc.maxm = 100

    E=fc.gs_energy(mode=mode)
    print("GS energy with ",mode,E)

    # compute particle number
    N=0
    for i in range(n):
        j=i
        v=fc.vev(fc.Cdagup[i]*fc.Cup[j]).real
        N=N+v
        v=fc.vev(fc.Cdagdn[i]*fc.Cdn[j]).real
        N=N+v
    print("N=",N)


    # compute the 1RDM
    #FIXME hermitian
    D=np.zeros((2*n,2*n),dtype=np.complex_)
    for i in range(n):
        for si in range(2):
            for j in range(n):
                for sj in range(2):
                    v=0
                    if si==0 and sj==0:
                        v=fc.vev(fc.Cdagup[i]*fc.Cup[j])
                    if si==0 and sj==1:
                        v=fc.vev(fc.Cdagup[i]*fc.Cdn[j])
                    if si==1 and sj==0:
                        v=fc.vev(fc.Cdagdn[i]*fc.Cup[j])
                    if si==1 and sj==1:
                        v=fc.vev(fc.Cdagdn[i]*fc.Cdn[j])
                    D[2*i+si,2*j+sj]=v

    W=np.zeros(n,dtype=np.complex_)
    # compute the double occupancy
    for i in range(n):
        v=fc.vev(fc.Nup[i]*fc.Ndn[i])
        W[i]=U*v

    return [E,D,W]





def hubbard_2d(n,t,U,mu,Lint,mode="ED"):
    """
    Constructs and solves a 2D Hubbard system with
    periodic boundary conditions using ED or DMRG.

    Inputs:
    `n`: total # of sites
    `t`: nn-hopping amplitude
    `U`: interaction amplitude
    `mu`: chemical potential
    `Lint`: first `Lint` sites will host a Hubbard interaction
    `mode` (optional): either "ED" or "DMRG"

    Output:
    `[E,D,W]` where
    `E`: groundstate energy
    `D`: one-particle reduced density matrix
    `W`: double occupancy
    """
    L=int(math.sqrt(n))
    fc = fermionchain.Spinful_Fermionic_Chain(n)

    # periodic first neighbor hopping
    h = 0
    for a in range(L):
      for b in range(L):
        i=a*L+b
        j=((a+1)%L)*L+(b+0)%L
        print("hopping",a,b,i,j)
        h = h + t*fc.Cdagup[i]*fc.Cup[j]
        h = h + t*fc.Cdagdn[i]*fc.Cdn[j]
        j=((a+0)%L)*L+(b+1)%L
        print("hopping",a,b,i,j)
        h = h + t*fc.Cdagup[i]*fc.Cup[j]
        h = h + t*fc.Cdagdn[i]*fc.Cdn[j]
    h = h + h.get_dagger() # Make Hermitian

    # Hubbard term
    for i in range(Lint):
      h = h + U*fc.Nup[i]*fc.Ndn[i]

    # chemical potential
    for i in range(n):
      h =h + mu*(fc.Nup[i]+fc.Ndn[i])

    fc.set_hamiltonian(h) # initialize the Hamiltonian
    print("WARNING: the chosen maximal bond dimension is most likely to small for a 2d-Hubbard model")
    fc.maxm = 100

    E=fc.gs_energy(mode=mode)
    print("GS energy with ",mode,E)

    # compute particle number
    N=0
    for i in range(n):
        j=i
        v=fc.vev(fc.Cdagup[i]*fc.Cup[j]).real
        N=N+v
        v=fc.vev(fc.Cdagdn[i]*fc.Cdn[j]).real
        N=N+v
    print("N=",N)


    # compute the 1RDM
    #FIXME hermitian
    D=np.zeros((2*n,2*n),dtype=np.complex_)
    for i in range(n):
        for si in range(2):
            for j in range(n):
                for sj in range(2):
                    v=0
                    if si==0 and sj==0:
                        v=fc.vev(fc.Cdagup[i]*fc.Cup[j])
                    if si==0 and sj==1:
                        v=fc.vev(fc.Cdagup[i]*fc.Cdn[j])
                    if si==1 and sj==0:
                        v=fc.vev(fc.Cdagdn[i]*fc.Cup[j])
                    if si==1 and sj==1:
                        v=fc.vev(fc.Cdagdn[i]*fc.Cdn[j])
                    D[2*i+si,2*j+sj]=v

    W=np.zeros(n,dtype=np.complex_)
    # compute the double occupancy
    for i in range(n):
        v=fc.vev(fc.Nup[i]*fc.Ndn[i])
        W[i]=U*v

    return [E,D,W]
