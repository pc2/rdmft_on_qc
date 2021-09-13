import os
import sys
import networkx as nx
import matplotlib
import copy
import matplotlib.pyplot as plt
from qiskit.opflow.primitive_ops import PauliOp
from qiskit.opflow.state_fns import CircuitStateFn
from qiskit.quantum_info import Pauli
from networkx.algorithms import approximation

def disjoint_qubits(a,b):
    disjoint=True
    for i in range(len(a)):
        if a[i]!="I" and b[i]!="I":
            return False
    return True

def disjoint_string(a,b):
    disjoint=True
    for i in range(len(a)):
        if not (a[i]==b[i] or (a[i]=="I" and b[i]!="I") or (a[i]!="I" and b[i]=="I")):
            return False
    return True

def commutecheck(commutemode,a,b):
  if commutemode=="disjointqubits":
    if disjoint_qubits(a,b):
      return True
  elif commutemode=="qubitwise":
    if disjoint_string(a,b):
      return True
  elif commutemode=="commute":
    if Pauli(a).commutes(b):
      return True
  elif commutemode=="anticommute":
    if Pauli(a).anticommutes(b):
      return True
  else:
    raise RuntimeError("ERROR: unknown commutemode")
  return False

def cliquecover(elements,commutemode="qubitwise",plotgraph=False,printcliques=False):
    #takes strings of Paulis in elements and returns a minimal set of cliques
    print("graphclique-cliquecover:",len(elements),"with commutationmode",commutemode)

    #build networkx graph
    g = nx.Graph()
    g.add_nodes_from(range(len(elements)))

    for i in range(len(elements)):
        for j in range(len(elements)):
          if commutecheck(commutemode,elements[i],elements[j]):
            g.add_edge(i,j)

    if plotgraph:
      f = plt.figure()
      nx.draw(g, ax=f.add_subplot(111),with_labels=True)
      #nx.drawing.nx_pydot.write_dot(g,"graph.dot")
      f.savefig("graph_"+commutemode+".png")                
    
    #the clique covering can be obtained from the coloring of the complement graph
    #(clique covering number=chromatic number of complement graph)

    g2=nx.complement(g)
    nx.drawing.nx_pydot.write_dot(g2,"graph2.dot")
    if plotgraph:
      f = plt.figure()
      nx.draw(g2, ax=f.add_subplot(111),with_labels=True)
      f.savefig("graph_"+commutemode+"_complement.png")                

    #graph coloring
    lmin=1000000000
    cmin=[]
    cliques2=[]
    for strat in ["largest_first","random_sequential","smallest_last","independent_set","connected_sequential_bfs","connected_sequential_dfs","saturation_largest_first"]:
      c=nx.greedy_color(g2,strategy=strat)
      cliques=[]
      for i in range(len(elements)):
          cliques.append([])
      for i in range(len(elements)):
        cliques[c[i]].append(elements[i])
      cliques2=[]
      for i in range(len(elements)):
        if len(cliques[i])>0:
          cliques2.append(cliques[i])
      print("strategy",strat,len(cliques2))
      if len(cliques2)<lmin:
        cmin=copy.deepcopy(cliques2)
        lmin=len(cliques2)


    a=approximation.clique_removal(g)[1]
    print("clique_removal",len(a))

    #try to combine cliques
    comb=[]
    for i in range(len(cmin)):
      for j in range(i+1,len(cmin)):
        c=True
        for i1 in range(len(cmin[i])):
          for j1 in range(len(cmin[j])):
            if not commutecheck(commutemode,cmin[i][i1],cmin[j][j1]):
              c=False
              break
          if not c:
            break
        if c:
          print("combine clique",i,j)
    #cliques:
    print("number of cliques",len(cmin))
    if printcliques:
      for i in range(len(cmin)):
          print("clique",i,cmin[i])
    return cmin


