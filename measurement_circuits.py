import os
import copy
import math
import sys
from scipy import sparse
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, assemble, transpile,Aer,execute
from qiskit.quantum_info import Statevector,random_statevector
from qiskit.quantum_info import Pauli
import numpy as np
import networkx as nx
from qiskit import IBMQ, assemble, transpile,Aer
import itertools
import galois
from sympy.combinatorics import Permutation
from qiskit.opflow.primitive_ops import PauliOp
#import term_grouping
#from term_grouping import QWCCommutativity,FullCommutativity,genMeasureCircuit,NetworkX_approximate_clique_cover,BronKerbosch,BronKerbosch_pivot
#from generate_measurement_circuit import MeasurementCircuit,_get_measurement_circuit

def clique2stab(nq,c):
  #build stabilizer matrix for each clique
  stab=np.zeros((2*nq+1,nq), dtype=int)
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
  signs=[1]
  return {"mqc":mqc,"mqubits":[mqubit],"signs":signs}

def rk_gf2(mat):
  return np.linalg.matrix_rank(mat_gf2(mat))

def mat_gf2(mat):
  GF2 = galois.GF(2)
  mat2=GF2.Zeros((mat.shape[0],mat.shape[1]))
  for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
      mat2[i,j]=mat[i,j]
  return mat2

def vec_gf2(mat):
  GF2 = galois.GF(2)
  mat2=GF2.Zeros((mat.shape[0]))
  for i in range(mat.shape[0]):
    mat2[i]=mat[i]
  return mat2

def stab_H(nq,stab,ih):
  stab2=copy.deepcopy(stab)
  #swap x and z
  stab2[ih+nq,:]=stab[ih,:].copy()
  stab2[ih,:]=stab[ih+nq,:].copy()
  #r_i=r_i+x_i*z_i
  stab2[2*nq,:]=stab[2*nq,:].copy()+stab[ih,:].copy()*stab[ih+nq,:].copy()
  return stab2

def stab_s(nq,stab,c):
  stab2=copy.deepcopy(stab)
  stab2[2*nq,:]=stab[2*nq,:].copy()+stab[c,:].copy()*stab[c+nq,:].copy()
  #stab2[c,c]=0
  stab2[c,:]=stab[c,:].copy()+stab[nq+c,:].copy()
  return stab2

def stab_cnot(nq,stab,c,t):
  stab2=copy.deepcopy(stab)
  ones=vec_gf2(np.ones(nq, dtype=int))
  stab2[2*nq,:]=stab[2*nq,:].copy()+stab[c+nq,:].copy()*stab[t,:].copy()*(stab[t+nq,:].copy()+stab[c,:].copy()+ones)

  stab2[c,:]=stab[c,:].copy()+stab[t,:].copy()
  stab2[t+nq,:]=stab[t+nq,:].copy()+stab[c+nq,:].copy()
  return stab2

def stab_cz(nq,stab,c,t):
  stab2=copy.deepcopy(stab)
  stab2=stab_H(nq,stab2,t)
  stab2=stab_cnot(nq,stab2,c,t)
  stab2=stab_H(nq,stab2,t)
  #stab2[c,t]=0
  #stab2[t,c]=0
  return stab2

def stab_swap(nq,stab,i,j):
  stab2=copy.deepcopy(stab)
  stab2=stab_cnot(nq,stab2,i,j)
  stab2=stab_cnot(nq,stab2,j,i)
  stab2=stab_cnot(nq,stab2,i,j)
  #stab2[i,:]=stab[j,:].copy()
  #stab2[i+nq,:]=stab[j+nq,:].copy()
  #stab2[j,:]=stab[i,:].copy()
  #stab2[j+nq,:]=stab[i+nq,:].copy()
  return stab2

def inv_perm(p):
  ip=copy.deepcopy(list(p))
  for i in range(len(p)):
    ip[p[i]]=i
  return ip

def build_measurement_circuits_commute(nq,cc_in,config):
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
  tprint=config.getboolean('QC','print_construction_of_measurement_circuits')

  #first reduce to Pauli-z 
  mqubit=-1
  mqc=[]

  if tprint:
    print("input",cc_in)
  qreg = QuantumRegister(nq)
  creg = ClassicalRegister(nq)
  tqc=QuantumCircuit(qreg,creg)
  rand_sv=random_statevector(2**nq,seed=2344)

  #cc_in=['IIIYYX', 'XXYIXY', 'IIIZXI', 'YIYIXZ', 'XYYXZI', 'XZXIII']

  cc_reordered=[]
  for c in cc_in:  
      cc2=list(c)
      cc2.reverse()
      cc_reordered.append("".join(cc2))
  
  cc=copy.deepcopy(cc_reordered)
  if tprint:
    print("cc with qiskit ordering",cc)

  order=list(range(nq))
  orders=[list(range(nq))]
  if config.getboolean('QC','reorder_measurement_qubits'):
    orders=itertools.permutations(order)

  complexity_min=10000000
  mmin={}

  #orders=[(0, 1, 2, 4, 5, 3)]
  for o in orders:
    if tprint:
      print("try order: ",o)
    io=o[:] #inv_perm(o)

    #permute measurements with order o
    c=[]
    for ci in cc:
      c2=[]
      for i in range(nq):
        c2.append(ci[o[i]])
      c.append("".join(c2))
    if tprint:
      print("after qubit-reordering")
      print(c)

    stab=clique2stab(nq,c)
    if tprint:
      print(stab)
    if config.getboolean('QC','reorder_measurement_qubits_restricted'):
      exclude=False
      for i in range(len(c)):
        if not(stab[i,i]==1 or stab[nq+i,i]==1): 
          exclude=True
          break

      if exclude:
        if tprint:
          print("excluded because reorder_measurement_qubits_restricted=T and not every measurement qubit has a pauli-operator")
        continue

    rk1=rk_gf2(stab[0:2*nq,0:nq])
    if tprint:
      print("GF2-rank stab=",rk1)
    rk2=rk_gf2(stab[nq:2*nq,0:nq])
    if tprint:
      print("GF2-rank X=",rk2)

    initialH=[]
    initialH_found=False
    #make X-matrix maximal rank by adding Hs
    for p in itertools.product([0,1],repeat=nq):
      if initialH_found and not config.getboolean('QC','try_all_initial_H_combinations'):
        break
      mqc=QuantumCircuit(qreg,creg)
      tqc=QuantumCircuit(qreg,creg)
      stab2=clique2stab(nq,c)
      stabgf2=mat_gf2(stab2)

      tqc.initialize(rand_sv.data,list(range(nq)))

      if tprint:
        print("initial")
        print(stabgf2)
      for i in range(nq):
        if p[i]==1:
          stabgf2=stab_H(nq,stabgf2,i)

      #apply 
      rk1H=rk_gf2(stabgf2[0:2*nq,0:nq])
      rk2H=rk_gf2(stabgf2[nq:2*nq,0:nq])
      if tprint:
        print(p,"GF2-ranks=",rk1H,rk2H)
      if rk2H==rk1:
        initialH=p[:]
        initialH_found=True
      else:
        continue

      print("maximal X-rank with Hs at ",initialH)
      if tprint:
        print(stabgf2)
        
      mqc=QuantumCircuit(qreg,creg)
      for i in range(nq):
        if initialH[i]==1:
          mqc.h(io[i])
      
      #reduce X-matrix to unit matrix
      GF2 = galois.GF(2)
      L,U,P=GF2.lup_decompose(stabgf2[nq:2*nq,:])
      if tprint:
        print(np.array_equal(P @ stabgf2[nq:2*nq,:], L @ U))
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
      if tprint:
        print(perm)

      p = Permutation(perm)
      ts = p.transpositions()
      if tprint:
        print(ts)

      for t in ts:
        if tprint:
          print("swap(",t[0],",",t[1],")")
        mqc.swap(io[t[0]],io[t[1]])
        stabgf2=stab_swap(nq,stabgf2,t[0],t[1]).copy()

      if tprint:
        print("after permutation")
        print(stabgf2)

      Linv=np.linalg.inv(L)
      for i in range(nq-1,-1,-1):
        for j in range(nq):
          if i!=j:
            if Linv[j,i]==1:
              if tprint:
                print("cnot(",i,",",j,")")
              mqc.cnot(io[i],io[j])
              stabgf2=stab_cnot(nq,stabgf2,i,j).copy()
      if tprint:
        print(stabgf2)
      
      #check if X-matrix is diagonal
      isdiag=True
      isU=True
      for i in range(nq):
        for j in range(nq):
          if i!=j:
            if stabgf2[i+nq,j]!=0:
              isdiag=False
          if j<i:
            if stabgf2[i+nq,j]!=0:
              isU=False
      if tprint:
        print("X is upper triangular :",isU)
      if not isU:
        raise RuntimeError("X-matrix is not upper triangular. Something went wrong.")
      
      if tprint:
        print("X is diag:",isdiag)
      if not isdiag:
        #reduce X to diagonal form
        for i in range(nq-2,-1,-1):
          if stabgf2[i+nq,i]!=0:
            if tprint:
              print("starting at row ",i)
            for j in range(nq-1,i,-1):
              if stabgf2[i+nq,j]==1:
                if tprint:
                  print("cnot(",j,",",i,")")
                mqc.cnot(io[j],io[i])
                stabgf2=stab_cnot(nq,stabgf2,j,i).copy()
      if tprint:
        print("final after reduction")
        print(stabgf2)

      #reduce Z to zero matrix with CZ and S
      for i in range(nq):
        for j in range(nq):
          if stabgf2[i,j]==1 and i!=j:
            if tprint:
              print("CZ(",i,",",j,")")
            mqc.cz(io[i],io[j])
            stabgf2=stab_cz(nq,stabgf2,i,j).copy()
      for i in range(nq):
        if stabgf2[i,i]==1:
          if tprint:
            print("S(",i,")")
          mqc.s(io[i])
          stabgf2=stab_s(nq,stabgf2,i)

      if tprint:
        print(stabgf2)
        print("final after Hs")

      for i in range(nq):
        mqc.h(io[i])
        stabgf2=stab_H(nq,stabgf2,i).copy()

      if tprint:
        print(stabgf2)
      
      #get sign from phase row
      signs=np.ones(nq,dtype=int)
      for i in range(nq):
        if stabgf2[2*nq,i]==1:
          if config.getboolean("QC","add_Ys_instead_of_separate_signs"):
            mqc.y(io[i])
          else:
            signs[i]=-1

      for ic in range(len(cc)):
        if tprint:
          print(cc[ic],"is measured at",io[ic])
        mqc.measure(ic,ic)
      if tprint:
        print("depth=",mqc.depth())

      if tprint:
        print("signs=",signs)

      ##construction of https://arxiv.org/pdf/1907.13623.pdf
      #stab=clique2stab(nq,c)
      #q=_get_measurement_circuit(stab,nq)
      #tqc=tqc.compose(q.circuit)
      
      if config.getboolean('QC','test_construction_of_measurement_circuits'):
        tqc=tqc.compose(mqc)
        #test with example vector 
        backend = Aer.get_backend('aer_simulator_statevector')
        shots=int(config['QC']['shots'])
        seed=int(config['rnd']['seed'])
        jobs=execute(tqc,backend=backend,backend_properties=backend.properties(),shots=shots,seed_simulator=seed,seed_transpiler=seed)#,optimization_level=3)
        jobs.wait_for_final_state(wait=0.05)

        if jobs.done():
          res=jobs.result().results
          if tprint:
            print("result=",res[0].data.counts)

        for r in res[0].data.counts:
          if tprint:
            print(r,bin(int(r, base=16))[2:].zfill(nq),res[0].data.counts[r]/shots)

        for ic in range(len(cc)):
          if tprint:
            print(cc_in[ic],"is measured as Z at",io[ic])
          V=0
          for r in res[0].data.counts:
            b=bin(int(r, base=16))[2:].zfill(nq)
            v=res[0].data.counts[r]/shots
            if b[nq-1-io[ic]]=='1':
              V=V-v
            else:
              V=V+v
          V=V*signs[ic]
          
          #check result
          mop=PauliOp(Pauli(cc_in[ic]))
          msparse=sparse.csr_matrix(mop.to_matrix())
          exp=np.dot(np.conj(rand_sv.data),msparse.dot(rand_sv.data)).real
          print("value=",cc_in[ic],V,exp,abs(V-exp))
      mqubits=[]
      for ic in range(len(cc)):
        mqubits.append(io[ic])
      
      m={"mqc":mqc,"mqubits":mqubits,"signs":signs}

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
      print("ordering",c,"constructed_complexity=",constructed_complexity,"transpiled_complexity=",transpiled_complexity," (complexity=",complexity_measure,")")
      if complexity<complexity_min:
          complexity_min=complexity
          #mqubit=mq
          mmin=copy.deepcopy(m)

  return mmin


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

  return m

