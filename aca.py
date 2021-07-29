import os
import copy
import math
import sys
import numpy as np
import matplotlib.pyplot as plt

def printmat(n,m,name,A):
    f=10
    B=np.zeros((n*f,2*m*f),dtype=np.float)

    for i in range(n):
        for j in range(m):
            for k in range(f):
                for l in range(f):
                    B[i*f+k,j*f+l]=A.real[i,j]
                    B[i*f+k,m*f+j*f+l]=A.imag[i,j]

    amax=1.0 #np.amax(B)
    B=abs(B)/(amax+1e-6)
    array = np.reshape(B, (f*n, 2*f*m))
    plt.imsave(name+".png", array,cmap="Greys")
    printmat_real(n,m,name,A)
    printmat_real_spin(n,m,name,A)

def printmat_real(n,m,name,A):
    f=10
    B=np.zeros((n*f,m*f),dtype=np.float)

    for i in range(n):
        for j in range(m):
            for k in range(f):
                for l in range(f):
                    B[i*f+k,j*f+l]=A.real[i,j]

    amax=1.0 #np.amax(B)
    B=abs(B)/(amax+1e-6)
    array = np.reshape(B, (f*n, f*m))
    plt.imsave(name+"_real.png", array,cmap="Greys")

def printmat_real_spin(n,m,name,A):
    f=10
    B=np.zeros((int(n*f/2),int(m*f/2)),dtype=np.float)

    for i in range(int(n/2)):
        for j in range(int(m/2)):
            for k in range(f):
                for l in range(f):
                    B[i*f+k,j*f+l]=A.real[2*i,2*j]

    amax=1.0 #np.amax(B)
    B=abs(B)/(amax+1e-6)
    array = np.reshape(B, (int(f*n/2), int(f*m/2)))
    plt.imsave(name+"_real_dn.png", array,cmap="Greys")
    
    f=10
    B=np.zeros((int(n*f/2),int(m*f/2)),dtype=np.float)

    for i in range(int(n/2)):
        for j in range(int(m/2)):
            for k in range(f):
                for l in range(f):
                    B[i*f+k,j*f+l]=A.real[2*i+1,2*j+1]

    amax=1.0 #np.amax(B)
    B=abs(B)/(amax+1e-6)
    array = np.reshape(B, (int(f*n/2), int(f*m/2)))
    plt.imsave(name+"_real_up.png", array,cmap="Greys")

def aca_reorder(norb,orbinteract,Din,Win):
    Dout=np.zeros((norb,norb),dtype=np.complex_)
    mapping=[]

    I=0
    for i in range(len(orbinteract)):
        mapping.append(orbinteract[i])

    for i in range(norb):
        found=False
        for j in mapping:
            if j==i:
                found=True
                break
        if not found:
            mapping.append(i)
    print("mapping:",mapping)
    for i in range(norb):
        for j in range(norb):
            Dout[i,j]=Din[mapping[i],mapping[j]]
    Wout=[]
    for i in Win:
        Wout.append(mapping[i])
    return [Dout,Wout]

def aca(norb,ninteract,Din):
    Dout=np.zeros((norb,norb),dtype=np.complex_)
    #print("norb=",norb,np.shape(Din)) 
    #print("ninteract=",ninteract) 
    l=1
    E=np.zeros((ninteract,norb-ninteract),dtype=np.complex_)
    E[0:ninteract,0:norb-ninteract]=Din[0:ninteract,ninteract:norb]
    #printmat(ninteract,norb-ninteract,"ACA_E_",E)
    #printmat(norb,norb,"ACA_D_",Din)

    [q,r] = np.linalg.qr(np.transpose(E), mode='complete')
    #printmat(norb-ninteract,norb-ninteract,"ACA_Q_",q)
    #printmat(norb-ninteract,ninteract,"ACA_R_",r)
    
    U=np.zeros((norb,norb),dtype=np.complex_)
    for i in range(ninteract):
        U[i,i]=1.0
    U[ninteract::,ninteract::]=q
    Dout=np.matmul(np.matmul(np.transpose(U),Din),U)
    Iout=np.matmul(np.transpose(U),U)

    for i in range(norb):
        Iout[i,i]=Iout[i,i]-1
    #print("Unitary check=",np.max(np.absolute(Iout)))
    #printmat(norb,norb,"ACA_I_",Iout)
    #printmat(norb,norb,"ACA_D2_",Dout)     
    return Dout

