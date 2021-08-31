import os
import copy
import math
import sys
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
import numpy as np
import networkx as nx
from qiskit import IBMQ, assemble, transpile,Aer
import term_grouping
from term_grouping import QWCCommutativity,FullCommutativity,genMeasureCircuit,NetworkX_approximate_clique_cover,BronKerbosch,BronKerbosch_pivot
from generate_measurement_circuit import MeasurementCircuit,_get_measurement_circuit
import itertools
import galois
from sympy.combinatorics import Permutation

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
        stab[nq+j,i]=1
      if t[j]=="Y":
        stab[j,i]=1
        stab[nq+j,i]=1
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
    elif mode=="gate_weights":
        ngq=0
        for gw in gatew:
            try:
                ngq+=gatew[gw]*qc.count_ops()[gw]
            except KeyError:
                ngq+=0
        return ngq
    else:
      raise RuntimeError("complexity_measure "+str(mode)+" not known")

def build_measurement_circuits_none(nq,cc,config):
  #options
  transpiler_gates=config["QC"]["transpiler_gates"].split(",")
  transpiler_seed=int(config["QC"]["transpiler_seed"])
  transpiler_couplings=[]
  for c in config["QC"]["transpiler_couplings"].split(","):
    transpiler_couplings.append([int(c.split("_")[0]),int(c.split("_")[1])])
  gatew=dict()
  for c in config["QC"]["gate_weights"].split(","):
    gatew[c.split("_")[0]]=int(c.split("_")[1])
  complexity_measure=config["QC"]["complexity_measure"]
  criterion_for_qc_optimality=config["QC"]["criterion_for_qc_optimality"]

  #first reduce to Pauli-z 
  cc2,preqc=paulis_to_zs(nq,cc)
  mqubit=-1
  mqc=[]
  tqc=[]

  #then reduce Pauli-z-string to single-z at some qubit if necessary
  if cc2[0].count("Z")>1:
      print("reduce ",cc2[0],"to single qubit")
      print("possibilities according to heuristic:")
      #find cnot gates so that number of gates is minimal
      complexity_min=float('inf')

      #reduce z's with cnots that are available in the coupling topology
      #build up reduction graph
      reduction_graph = nx.DiGraph()
      root="".join(cc2[0])
      reduction_graph.add_node(root)
      reductions=[]

      single_found=False

      while not single_found:
          reduction_found=True
          while reduction_found:
              reduction_found=False
              #iterate over leaves and try a local reduction for every leaf
              for cl in [v for v, d in reduction_graph.out_degree() if d == 0]:
                  for coupling in transpiler_couplings:
                      if cl[coupling[0]]=="Z" and cl[coupling[1]]=="Z":
                          #two reduction directions are possible
                          c=list(cl[:])
                          c[coupling[0]]="I"
                          if "".join(c) not in [v for v, d in reduction_graph.out_degree() if d == 0]:
                              #FIXME check direction
                              reduction_graph.add_node("".join(c))
                              reduction_graph.add_edge(cl,"".join(c),weight=1,reduction={"op":"cnot","c":coupling[0],"t":coupling[1]})
                              reduction_found=True
                          c=list(cl[:])
                          c[coupling[1]]="I"
                          if "".join(c) not in [v for v, d in reduction_graph.out_degree() if d == 0]:
                              #FIXME check direction
                              reduction_graph.add_node("".join(c))
                              reduction_graph.add_edge(cl,"".join(c),weight=1,reduction={"op":"cnot","c":coupling[1],"t":coupling[0]})
                              reduction_found=True

          leaves= [node for node in reduction_graph.nodes() if reduction_graph.in_degree(node)!=0 and reduction_graph.out_degree(node)==0]                    
          #check if a final reduction was already found
          for l in leaves:
              if l.count('Z')==1:
                  single_found=True
                  exit

          if not single_found:
              #try swaps that are composed of cnots that are available in the coupling
              #iterate over leaves
              for cl in [v for v, d in reduction_graph.out_degree() if d == 0]:
                  for coupling in transpiler_couplings:
                      if (cl[coupling[0]]=="I" and cl[coupling[1]]=="Z") or (cl[coupling[0]]=="Z" and cl[coupling[1]]=="I"):
                          c=list(cl[:])
                          tmp=c[coupling[0]]
                          c[coupling[0]]=c[coupling[1]]
                          c[coupling[1]]=tmp
                          if "".join(c) not in [v for v, d in reduction_graph.out_degree() if d == 0]:
                              reduction_graph.add_node("".join(c))
                              reduction_graph.add_edge(cl,"".join(c),weight=3,reduction={"op":"swap","c":coupling[0],"t":coupling[1]})

          
          leaves= [node for node in reduction_graph.nodes() if reduction_graph.in_degree(node)!=0 and reduction_graph.out_degree(node)==0]                    
          #check if a final reduction was already found
          for l in leaves:
              if l.count('Z')==1:
                  single_found=True
                  exit

      leaves= [node for node in reduction_graph.nodes() if reduction_graph.in_degree(node)!=0 and reduction_graph.out_degree(node)==0]                    
      #check if a final reduction was already found
      for l in leaves:
          if l.count('Z')==1:
              q = QuantumRegister(nq)
              c = ClassicalRegister(nq)
              qc=QuantumCircuit(q,c)
              qc=qc.compose(preqc)
              pa=nx.dijkstra_path(reduction_graph,root,l,weight='weight')
              for i in range(len(pa)-1):
                  d=reduction_graph.get_edge_data(pa[i],pa[i+1])
                  if d["reduction"]["op"]=="cnot":
                      qc.cnot(d["reduction"]["c"],d["reduction"]["t"])
                  if d["reduction"]["op"]=="swap":
                      qc.swap(d["reduction"]["c"],d["reduction"]["t"])
              mq=l.index('Z')
              qc.measure(mq,0)
                      
              transpiled_qc = transpile(qc, basis_gates=transpiler_gates,coupling_map=transpiler_couplings, optimization_level=3,seed_transpiler=transpiler_seed)
              
              constructed_complexity=measure_complexity(qc,mode=complexity_measure,gate_weights=gatew)
              transpiled_complexity=measure_complexity(transpiled_qc,mode=complexity_measure,gate_weights=gatew)
              complexity=0
              if criterion_for_qc_optimality=="constructed":
                  complexity=constructed_complexity
              elif criterion_for_qc_optimality=="transpiled":
                  complexity=transpiled_complexity
              else:
                  raise RuntimeError('criterion_for_qc_optimality not known')
              print("root=",root,"target=",l,"constructed_complexity=",constructed_complexity,"transpiled_complexity=",transpiled_complexity," (complexity=",complexity_measure,")")
              #print(transpiled_qc)
              if complexity<complexity_min:
                  complexity_min=complexity
                  mqubit=mq
                  mqc=copy.deepcopy(qc)
  else:
      mqubit=cc2[0].index("Z")
      q = QuantumRegister(nq)
      c = ClassicalRegister(nq)
      mqc=QuantumCircuit(q,c)
      mqc=mqc.compose(preqc)
      mqc.measure(mqubit,0)
  return {"mqc":mqc,"mqubits":[mqubit]}

def rk_gf2(mat):
  return np.linalg.matrix_rank(mat_gf2(mat))

def mat_gf2(mat):
  GF2 = galois.GF(2)
  mat2=GF2.Zeros((mat.shape[0],mat.shape[1]))
  for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
      mat2[i,j]=mat[i,j]
  return mat2

def stab_H(nq,stab,ih):
  stab2=copy.deepcopy(stab)
  stab2[ih+nq,:]=stab[ih,:].copy()
  stab2[ih,:]=stab[ih+nq,:].copy()
  return stab2

def stab_cnot(nq,stab,c,t):
  stab2=copy.deepcopy(stab)
  stab2[c,:]=stab[c,:].copy()+stab[t,:].copy()
  stab2[c+nq,:]=stab[c+nq,:].copy()+stab[t+nq,:].copy()
  return stab2

def stab_swap(nq,stab,i,j):
  stab2=copy.deepcopy(stab)
  stab2[i,:]=stab[j,:].copy()
  stab2[i+nq,:]=stab[j+nq,:].copy()
  stab2[j,:]=stab[i,:].copy()
  stab2[j+nq,:]=stab[i+nq,:].copy()
  return stab2

def build_measurement_circuits_commute(nq,cc,config):
  #options
  transpiler_gates=config["QC"]["transpiler_gates"].split(",")
  transpiler_seed=int(config["QC"]["transpiler_seed"])
  transpiler_couplings=[]
  for c in config["QC"]["transpiler_couplings"].split(","):
    transpiler_couplings.append([int(c.split("_")[0]),int(c.split("_")[1])])
  gatew=dict()
  for c in config["QC"]["gate_weights"].split(","):
    gatew[c.split("_")[0]]=int(c.split("_")[1])
  complexity_measure=config["QC"]["complexity_measure"]
  criterion_for_qc_optimality=config["QC"]["criterion_for_qc_optimality"]

  #first reduce to Pauli-z 
  mqubit=-1
  mqc=[]
  tqc=[]

  print(cc)
  #cc,preqc=paulis_to_zs(nq,cc)
  #variants=itertools.permutations(cc)
  #cc.reverse()


  variants=[cc]

  complexity_min=10000000
  mmin={}
  for c in variants:
    print("variant: ",c)
    stab=clique2stab(nq,c)
    print(stab)

    rk1=rk_gf2(stab)
    print("GF2-rank stab=",rk1)
    rk2=rk_gf2(stab[nq:])
    print("GF2-rank X=",rk2)

    initialH=[]
    #make X-matrix maximal rank by adding Hs
    print(itertools.product([0,1],repeat=nq))
    for p in itertools.product([0,1],repeat=nq):
      #print(p)
      stab2=clique2stab(nq,c)
      for i in range(nq):
        if p[i]==1:
          stab2=stab_H(nq,stab2,i)
      #print(stab2)

      #apply 
      rk1H=rk_gf2(stab2)
      rk2H=rk_gf2(stab2[nq:])
      #print("GF2-rank stab=",rk1H)
      #print("GF2-rank x=",rk2H)
      if rk2H==rk1:
        initialH=p[:]
        break
    print("maximal X-rank with Hs at ",initialH)
    print(stab2)
    
    #reduce X-matrix to unit matrix
    GF2 = galois.GF(2)
    L,U,P=GF2.lup_decompose(mat_gf2(stab2[nq:]))
    print(np.array_equal(P @ mat_gf2(stab2[nq:]), L @ U))
    print("L=")
    print(L)
    print("U=")
    print(U)
    print("P=")
    print(P)
    
    
    #permute columns
    perm=[]
    for i in range(nq):
      for j in range(nq):
        if P[i,j]==1:
          perm.append(j)
    print(perm)

    p = Permutation(perm)
    ts = p.transpositions()
    print(ts)

    stab3=copy.deepcopy(stab2)
    for t in ts:
      print("swap(",t[0],",",t[1],")")
      stab3=stab_swap(nq,stab3,t[0],t[1]).copy()

    print("after permutation")
    print(stab3)

    Linv=np.linalg.inv(L)
    print("Linv")
    print(Linv)
    print("Linv P A")
    print(Linv @ P @ mat_gf2(stab2[nq:]))
    print("Linv P A")
    print(Linv @ mat_gf2(stab3[nq:]))
    #check reduction
    stabgf2=mat_gf2(stab3)
    for i in range(nq-1,-1,-1):
      for j in range(nq):
        if i!=j:
          if Linv[i,j]==1:
            print("cnot(",i,",",j,")")
            stabgf2=stab_cnot(nq,stabgf2,i,j)
    print(stabgf2)
    
    #for i in range(nq):
    #  stabgf2=stab_H(nq,stabgf2,i)
    #print(stabgf2)





  return None 
  try:
      q=_get_measurement_circuit(stab,nq)
  except AssertionError:
      print("_get_measurement_circuit has failed")
  m={"mqc":q.circuit,"mqubits":[]}
  
  transpiled_mqc = transpile(m["mqc"], basis_gates=transpiler_gates,coupling_map=transpiler_couplings, optimization_level=3,seed_transpiler=transpiler_seed)
  constructed_complexity=measure_complexity(m["mqc"],mode=complexity_measure,gate_weights=gatew)
  transpiled_complexity=measure_complexity(transpiled_mqc,mode=complexity_measure,gate_weights=gatew)
  complexity=0
  if criterion_for_qc_optimality=="constructed":
      complexity=constructed_complexity
  elif criterion_for_qc_optimality=="transpiled":
      complexity=transpiled_complexity
  else:
      raise RuntimeError('criterion_for_qc_optimality not known')
  print("variant",c,"constructed_complexity=",constructed_complexity,"transpiled_complexity=",transpiled_complexity," (complexity=",complexity_measure,")")
  #print(transpiled_qc)
  if complexity<complexity_min:
      complexity_min=complexity
      #mqubit=mq
      mmin=copy.deepcopy(m)
  print(mmin["mqc"])

  return m


def build_measurement_circuit(mode,nq,cc,config):
  #options
  transpiler_gates=config["QC"]["transpiler_gates"].split(",")
  transpiler_seed=int(config["QC"]["transpiler_seed"])
  transpiler_couplings=[]
  for c in config["QC"]["transpiler_couplings"].split(","):
    transpiler_couplings.append([int(c.split("_")[0]),int(c.split("_")[1])])
  gatew=dict()
  for c in config["QC"]["gate_weights"].split(","):
    gatew[c.split("_")[0]]=int(c.split("_")[1])

  complexity_measure=config["QC"]["complexity_measure"]
  criterion_for_qc_optimality=config["QC"]["criterion_for_qc_optimality"]

  if mode=="none" or len(cc)==1:
    #one measurement per program
    m=build_measurement_circuits_none(nq,cc,config)
  elif mode=="disjointqubits" or mode=="qubitwise" or mode=="commute":
    m=build_measurement_circuits_commute(nq,cc,config)

  if m is not None:
    transpiled_mqc = transpile(m["mqc"], basis_gates=transpiler_gates,coupling_map=transpiler_couplings, optimization_level=3,seed_transpiler=transpiler_seed)
    constructed_complexity=measure_complexity(m["mqc"],mode=complexity_measure,gate_weights=gatew)
    transpiled_complexity=measure_complexity(transpiled_mqc,mode=complexity_measure,gate_weights=gatew)

    print("constructed measurement circuit has constructed complexity",constructed_complexity,"transpiled_complexity=",transpiled_complexity," (complexity=",complexity_measure,")")
    print("measure "+str(cc)+" at",m["mqubits"])
            
    if config.getboolean('QC','print_measurement_circuits'):
      print("constrcuted measurement program")
      print(m["mqc"])
      print("transpiled measurement program")
      print(transpiled_mqc)

  quit()
  return m

