import os
import copy
import math
import sys
import numpy as np
import openfermion as of
from numpy.random import rand
from scipy.optimize import minimize,BFGS

norb=0
nvar=0
fop=[]
cop=[]
cval=[]
uval=0

def op_norm(norb):
    op=of.FermionOperator(((0, 1), (norb-1, 0)), coefficient=0.0)
    op+=of.FermionOperator((), coefficient=1.0)
    return of.get_sparse_operator(op)

def op_rho(norb,a,b,treal):
    if treal:
        coeff=0.5
        op=of.FermionOperator(((0, 1), (norb-1, 0)), coefficient=0.0)
        op+=of.FermionOperator(((a, 1), (b, 0)), coefficient=coeff)
        op+=of.FermionOperator(((b, 1), (a, 0)), coefficient=coeff)
        return of.get_sparse_operator(op)
    else:
        coeff=0.5/1j
        op=of.FermionOperator(((0, 1), (norb-1, 0)), coefficient=0.0)
        op+=of.FermionOperator(((a, 1), (b, 0)), coefficient=coeff)
        op+=of.FermionOperator(((b, 1), (a, 0)), coefficient=-coeff)
        return of.get_sparse_operator(op)

def op_hubbard(norb,Wsites):
    coeff=1.0
    op=of.FermionOperator(((0, 1), (norb-1, 0)), coefficient=0.0)
    for i in Wsites:
        op+=Uval*of.FermionOperator(((2*i, 1), (2*i, 0), (2*i+1, 1),(2*i+1, 0)), coefficient=coeff)
    return of.get_sparse_operator(op)

def obj_val(x):
    return (np.dot(x,fop.dot(x))).real

def cons_val(x):
    y=np.zeros((len(cop)),dtype=np.float64)
    for i in range(len(cop)):
        y[i]=(np.dot(x,cop[i].dot(x))).real-cval[i]
    return y

def obj_der(x):
    return fop.dot(x)+fop.transpose().dot(x)

def cons_der(x):
    J=[]
    for i in range(len(cop)):
        J.append(cop[i].dot(x)+cop[i].transpose().dot(x))
    return J;

def F_hubbard(norbin,U,orbinteract,D,options):
    print("WARNING: This is a completely UNoptimized simple variant for the constrained minimization, for larger number of orbitals please use the optimized Fortan implementation.")
    #constrained minimization for density-matrix functional
    F=0

    global norb
    global nvar
    global cop
    global cval
    global fop
    global Uval
    norb=norbin
    Uval=U
    

    nvar=2**norb
    if nvar>=4096:
        exit()

    print("nvar=",nvar)

    #build sparse matrix representation for operators
    #interaction
    Wlocal=[]
    for i in range(int(len(orbinteract)/2)):
        Wlocal.append(i)
    fop=op_hubbard(norb,Wlocal)
    
    cop=[]
    cval=[]
    #norm
    cop.append(op_norm(norb))
    cval.append(1.0)

    #density matrix
    for a in range(norb):
        for b in range(a,norb):
            if a==b:
                cop.append(op_rho(norb,a,b,True))
                cval.append(D[a,b].real)
            else:
                cop.append(op_rho(norb,a,b,True))
                cval.append(D[a,b].real)
                #cop.append(op_rho(norb,a,b,False))
                #cval.append(D[a,b].imag)
    #print(cop)
    #print(cval)
    
    x0=rand(nvar)
    x0=x0/math.sqrt(np.sum(x0**2))

    eq_cons={'type': 'eq',
           'fun' : lambda x: cons_val(x),
           'jac' : lambda x: cons_der(x)}
    res=minimize(obj_val, x0, jac=obj_der, method='SLSQP', constraints=[eq_cons],tol=options['tol'],options={'maxiter':options['maxiter'],'verbose': 0,'disp': True},hess=BFGS())
    F=res.fun
    return F
