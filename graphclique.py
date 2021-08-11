import os
import sys
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from qiskit.opflow.primitive_ops import PauliOp
from qiskit.opflow.state_fns import CircuitStateFn
from qiskit.quantum_info import Pauli

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

def cliquecover(elements,commutemode="qubitwise",plotgraph=False,printcliques=False):
    #takes strings of Paulis in elements and returns a minimal set of cliques
    print("graphclique-cliquecover:",len(elements),"with commutationmode",commutemode)

    #build networkx graph
    g = nx.Graph()
    g.add_nodes_from(range(len(elements)))

    if commutemode=="disjointqubits":
      for i in range(len(elements)):
          for j in range(len(elements)):
              if disjoint_qubits(elements[i],elements[j]):
                  g.add_edge(i,j)
    elif commutemode=="qubitwise":
      for i in range(len(elements)):
          for j in range(len(elements)):
              if disjoint_string(elements[i],elements[j]):
                  g.add_edge(i,j)
    elif commutemode=="commute":
      #for i in range(len(elements)):
      #  elements[i]=elements[i].replace("X","Z")
      #  elements[i]=elements[i].replace("Y","Z")
      for i in range(len(elements)):
        for j in range(len(elements)):
          if Pauli(elements[i]).commutes(elements[j]):
            g.add_edge(i,j)
    elif commutemode=="anticommute":
      #for i in range(len(elements)):
      #  elements[i]=elements[i].replace("X","Z")
      #  elements[i]=elements[i].replace("Y","Z")
      for i in range(len(elements)):
        for j in range(len(elements)):
          if Pauli(elements[i]).anticommutes(elements[j]):
            g.add_edge(i,j)
      
    else:
      print("ERROR: unknown commutemode",commutemode)
      exit()


    if plotgraph:
      f = plt.figure()
      nx.draw(g, ax=f.add_subplot(111),with_labels=True)
      f.savefig("graph_"+commutemode+".png")                
    
    #the clique covering can be obtained from the coloring of the complement graph
    #(clique covering number=chromatic number of complement graph)

    g2=nx.complement(g)
    if plotgraph:
      f = plt.figure()
      nx.draw(g2, ax=f.add_subplot(111),with_labels=True)
      f.savefig("graph_"+commutemode+"_complement.png")                

    #graph coloring
    c=nx.greedy_color(g2,strategy="largest_first")
#    print("coloring result:")
    cliques=[]
    for i in range(len(elements)):
        cliques.append([])

    for i in range(len(elements)):
#        print(i,elements[i],"->",c[i])
        cliques[c[i]].append(elements[i])

    cliques2=[]
    for i in range(len(elements)):
      if len(cliques[i])>0:
        cliques2.append(cliques[i])

    #cliques:
    if printcliques:
      for i in range(len(cliques2)):
          print("clique",i,cliques2[i])
    return cliques2
