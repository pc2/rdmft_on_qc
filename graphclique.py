import os
import sys
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt

def disjoint_string(a,b):
    disjoint=True
    for i in range(len(a)):
        if a[i]!="I" and b[i]!="I":
            return False
    return True

def cliquecover(elements):
    #takes strings of zeros and ones in elements and returns a minimal set of cliques
    print("graphclique-cliquecover:",elements)

    #build networkx graph
    g = nx.Graph()
    g.add_nodes_from(range(len(elements)))
    for i in range(len(elements)):
        for j in range(len(elements)):
            if disjoint_string(elements[i],elements[j]):
                g.add_edge(i,j)
    f = plt.figure()
    nx.draw(g, ax=f.add_subplot(111),with_labels=True)
    f.savefig("graph.png")                
    
    #the clique covering can be obtained from the coloring of the complement graph
    #(clique covering number=chromatic number of complement graph)

    g2=nx.complement(g)
    f = plt.figure()
    nx.draw(g2, ax=f.add_subplot(111),with_labels=True)
    f.savefig("graph2.png")                

    #graph coloring
    c=nx.greedy_color(g2,strategy="largest_first")
    print("coloring result:")
    cliques=[]
    for i in range(len(elements)):
        cliques.append([])

    for i in range(len(elements)):
        print(i,elements[i],"->",c[i])
        cliques[c[i]].append(elements[i])

    #cliques:
    for i in range(len(elements)):
        if len(cliques[i])>0:
            print("clique",i,cliques[i])
    return cliques
