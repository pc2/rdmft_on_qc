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
ctype=[]

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

def obj_val_cplx(x):
    return (np.dot(x[0:int(nvar/2)],fop.dot(x[0:int(nvar/2)]))).real+(np.dot(x[int(nvar/2)::],fop.dot(x[int(nvar/2)::]))).real

def cons_val(x):
    y=np.zeros((len(cop)),dtype=np.float64)
    for i in range(len(cop)):
        y[i]=(np.dot(x,cop[i].dot(x))).real-cval[i]
    return y

def cons_val_cplx(x):
    y=np.zeros((len(cop)),dtype=np.float64)
    for i in range(len(cop)):
        if ctype[i]!="rdm_imag":
            y[i]=(np.dot(x[0:int(nvar/2)],cop[i].dot(x[0:int(nvar/2)]))).real+(np.dot(x[int(nvar/2)::],cop[i].dot(x[int(nvar/2)::]))).real-cval[i]
        else:
            y[i]=(1j*(np.dot(x[0:int(nvar/2)],cop[i].dot(x[int(nvar/2)::])))-1j*(np.dot(x[int(nvar/2)::],cop[i].dot(x[0:int(nvar/2)]))))-cval[i]
    return y

def obj_der(x):
    return fop.dot(x)+fop.transpose().dot(x)

def obj_der_cplx(x):
    y=np.zeros((nvar),dtype=np.float64)
    y[0:int(nvar/2)]=fop.dot(x[0:int(nvar/2)])+fop.transpose().dot(x[0:int(nvar/2)])
    y[int(nvar/2)::]=fop.dot(x[int(nvar/2)::])+fop.transpose().dot(x[int(nvar/2)::])
    return y

def cons_der(x):
    J=[]
    for i in range(len(cop)):
        J.append(cop[i].dot(x)+cop[i].transpose().dot(x))
    return J;

def cons_der_cplx(x):
    J=[]
    for i in range(len(cop)):
        if ctype[i]!="rdm_imag":
            y=np.zeros((nvar),dtype=np.float64)
            y[0:int(nvar/2)]=cop[i].dot(x[0:int(nvar/2)])+cop[i].transpose().dot(x[0:int(nvar/2)])
            y[int(nvar/2)::]=cop[i].dot(x[int(nvar/2)::])+cop[i].transpose().dot(x[int(nvar/2)::])
            J.append(y)
        else:
            y=np.zeros((nvar),dtype=np.float64)
            y[0:int(nvar/2)]=1j*(cop[i].dot(x[int(nvar/2)::])-cop[i].transpose().dot(x[int(nvar/2)::]))
            y[int(nvar/2)::]=1j*(cop[i].transpose().dot(x[0:int(nvar/2)])-cop[i].dot(x[0:int(nvar/2)]))
            J.append(y)
    return J;

def F_hubbard(norbin,U,orbinteract,D,options,tcplx):
    print("WARNING: This is a completely UNoptimized simple variant for the constrained minimization, for larger number of orbitals please use the optimized Fortan implementation.")
    #constrained minimization for density-matrix functional
    F=0

    global norb
    global nvar
    global cop
    global cval
    global fop
    global Uval
    global ctype
    norb=norbin
    Uval=U
    
    nvar=2**norb
    print("tcplx=",tcplx)
    if tcplx:
        nvar=nvar*2
    if nvar>=4096:
        exit()

    print("nvar=",nvar)

    #build sparse matrix representation for operators
    #interaction
    Wsites=[]
    for i in range(int(len(orbinteract)/2)):
        Wsites.append(i)
    print("Wsites=",Wsites)
    fop=op_hubbard(norb,Wsites)
    
    cop=[]
    cval=[]
    #norm
    cop.append(op_norm(norb))
    cval.append(1.0)
    ctype.append("norm")

    #density matrix
    for a in range(norb):
        for b in range(a,norb):
            if a==b:
                cop.append(op_rho(norb,b,a,True))
                cval.append(D[a,b].real)
                ctype.append("rdm_real")
            else:
                cop.append(op_rho(norb,b,a,True))
                cval.append(D[a,b].real)
                ctype.append("rdm_real")
                if tcplx:
                    cop.append(op_rho(norb,b,a,False))
                    cval.append(D[a,b].imag)
                    ctype.append("rdm_imag")
    #for i in range(len(cop)):
    #    print(i,"cval=",cval[i])
    #    print(i,"cop=",cop[i])
    #print("cval=",cval)

    x0=rand(nvar)
    x0=x0/math.sqrt(np.sum(x0**2))
    
    if False:
        J=cons_der_cplx(x0)
        c=cons_val_cplx(x0)
        dx=1e-6
        for i in range(nvar):
            dx0=np.copy(x0)
            dx0[i]=dx0[i]+dx
            dc=cons_val_cplx(dx0)
            for j in range(len(cop)):
                print(j,i,J[j][i],(dc[j]-c[j])/dx)


    res=[]
    if tcplx:
        eq_cons={'type': 'eq',
           'fun' : lambda x: cons_val_cplx(x),
           'jac' : lambda x: cons_der_cplx(x)}
        res=minimize(obj_val_cplx, x0, jac=obj_der_cplx, method='trust-constr', constraints=[eq_cons],tol=options['tol'],options={'maxiter':options['maxiter'],'verbose': 2,'disp': True},hess=BFGS())
    else:
        eq_cons={'type': 'eq',
           'fun' : lambda x: cons_val(x),
           'jac' : lambda x: cons_der(x)}
        res=minimize(obj_val, x0, jac=obj_der, method='trust-constr', constraints=[eq_cons],tol=options['tol'],options={'maxiter':options['maxiter'],'verbose': 2,'disp': True},hess=BFGS())
    F=res.fun
    return F
