import sys
import os
from os import listdir
from os.path import isfile, join

mypath="./"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

o = open("data.dat", "w")
for f in sorted(onlyfiles):
    if f.find("combination_of_observables")>=0 and f.find(".out")>=0:
        g=f.split(".")[0]
        s=(int(g.split("_")[3]))*2
        m=g.split("_")[4]
        l=g.split("_")[5]
        fi = open(f, "r")
        v0=0
        v1=0
        v2=0
        for x in fi:
            if x.find("unique pauliops")==0:
                v0=x.split()[2]
            if x.find("number of cliques")==0:
                v1=x.split()[3]
#            if x.find("groups with https")>=0:
#                v2=x.split()[4]
        fi.close()
        o.write(f+" "+str(s)+" "+m+" "+l+" "+str(v0)+" "+str(v1)+"\n ")
        #+str(v1)+" "+str(v2)+"\n")

o.close()
