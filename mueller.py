import os
import copy
import math
import sys
import numpy as np
import scipy as sp

def mueller_hubbard(norb,ninteract,U,D):
    #see https://arxiv.org/pdf/1509.01985.pdf
    #see https://ediss.uni-goettingen.de/bitstream/handle/11858/00-1735-0000-002E-E5C2-7/out.pdf?sequence=1 eq. 5.22
    FM=0
    FH=0
    FF=0
    #[f,v]=np.linalg.eig(D)
    #print("f=",f)

    sqrtD=sp.linalg.sqrtm(D)
    #print("sqrtD",sqrtD)

    elem=[]
    for i in range(int(ninteract/2)):
        a=i
        b=i+1
        elem.append([U,a,a,a,a])
        elem.append([U,b,b,b,b])
        elem.append([U,a,b,a,b])
        elem.append([U,b,a,b,a])
    for e in elem:
        a=e[1]
        b=e[2]
        d=e[3]
        c=e[4]
        FH=FH+0.5*e[0]*D[d,a]*D[c,b]
        FF=FF+0.5*e[0]*(-D[c,a]*D[d,b])
        FM=FM+0.5*e[0]*(-sqrtD[c,a]*sqrtD[d,b])
    #print("F_M=",FH+FM,FH,FH+FF)
    return np.real(FH+FM)

def mueller_hubbard_der(norb,ninteract,U,D):
    #numerical derivatives of the mueller functional with respect to the 1rdm
    dx=1e-5
    der=[]
    for i in range(norb):
        for j in range(norb):
            if i>j:
                continue
            if i==j:
                D2=copy.deepcopy(D)
                D2[i,j]+=dx
                Fp=mueller_hubbard(norb,ninteract,U,D2)
                D2=copy.deepcopy(D)
                D2[i,j]-=dx
                Fm=mueller_hubbard(norb,ninteract,U,D2)
                print("dFM_dD real",i,j,0.5*(Fp-Fm)/dx)
                der.append(0.5*(Fp-Fm)/dx)
            else:
                D2=copy.deepcopy(D)
                D2[i,j]+=dx
                D2[j,i]+=dx
                Fp=mueller_hubbard(norb,ninteract,U,D2)
                D2=copy.deepcopy(D)
                D2[i,j]-=dx
                D2[j,i]-=dx
                Fm=mueller_hubbard(norb,ninteract,U,D2)
                print("dFM_dD real",i,j,0.5*(Fp-Fm)/dx)
                der.append(0.5*(Fp-Fm)/dx)
            if i!=j:
                D2=copy.deepcopy(D)
                D2[i,j]+=dx*1j
                D2[j,i]-=dx*1j
                Fp=mueller_hubbard(norb,ninteract,U,D2)
                D2=copy.deepcopy(D)
                D2[i,j]-=dx*1j
                D2[j,i]+=dx*1j
                Fm=mueller_hubbard(norb,ninteract,U,D2)
                print("dFM_dD imag",i,j,0.5*(Fp-Fm)/dx)
                der.append(0.5*(Fp-Fm)/dx)
    return der


