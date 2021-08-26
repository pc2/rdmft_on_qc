import os
import copy
import math
import sys
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
import numpy as np

def clique2stab(nq,c):
  #build stabilizer matrix for each clique
  stab=np.zeros((2*nq,nq), dtype=int)
  i=-1
  for t in c:
    i=i+1
    if i>=nq:
      raise RuntimeError("clique"+str(c)+"has too many members")
    for j in range(len(t)):
      if t[j]=="Z":
        stab[j,i]=1
      if t[j]=="X":
        stab[norb_aca+j,i]=1
      if t[j]=="Y":
        stab[j,i]=1
        stab[norb_aca+j,i]=1
  return stab

def paulis_to_zs(nq,ops):
  q = QuantumRegister(nq)
  c = ClassicalRegister(nq)
  qc=QuantumCircuit(q,c)

  ops2=[]
  for op in ops:
    ops2.append(list(op))
  #map to only sigma_z
  sg=[]
  for i in range(nq):
    isX=False
    isY=False
    for j in range(len(ops2)):
      if ops2[j][i]=='X':
        if isY:
          print("error")
          quit()
        isX=True
      if ops2[j][i]=='Y':
        if isX:
          print("error")
          quit()
        isY=True
    if isX:
      qc.h(i)
    if isX:
      qc.sdg(i)
      qc.h(i)

    for j in range(len(ops2)):
      if ops2[j][i]=='X' and isX:
        ops2[j][i]='Z'
      if ops2[j][i]=='Y' and isY:
        ops2[j][i]='Z'

  return ops2,qc

def measure_complexity(qc,mode="depth",gate_weights=None):
    if mode=="depth":
        return qc.depth()
    else:
        ngq=0
        for gw in gatew:
            try:
                ngq+=gatew[gw]*qc.count_ops()[gw]
            except KeyError:
                ngq+=0
        return ngq
