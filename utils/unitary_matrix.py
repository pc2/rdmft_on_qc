import numpy as np
from qiskit import *
from qiskit import Aer
from qiskit.circuit import Gate
from math import pi
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Unroller

#Changing the simulator 
backend = Aer.get_backend('unitary_simulator')

#The circuit without measurement
A = Gate('A', 1, [])
B = Gate('B', 1, [])
C = Gate('C', 1, [])
alpha = 1
c=0
t=1

theta = 0.1 # theta can be anything (pi chosen arbitrarily)
circ = QuantumCircuit(2)

circ.sdg(0)
circ.h(0)
circ.sdg(1)
circ.h(1)
circ.cnot(1,0)
circ.draw()
print(circ)
print(circ.decompose())

#job execution and getting the result as an object
job = execute(circ, backend)
result = job.result()

pass_ = Unroller(['u1', 'u2', 'u3', 'cx'])
pm = PassManager(pass_)
new_circ = pm.run(circ)
print(new_circ)

#get the unitary matrix from the result object
u1=result.get_unitary(circ, decimals=3)

circ = QuantumCircuit(2)

circ.cnot(1,0)
circ.cz(0,1)
circ.sdg(1)
circ.h(0)
circ.h(1)
circ.draw()
print(circ)

#job execution and getting the result as an object
job = execute(circ, backend)
result = job.result()

#get the unitary matrix from the result object
u2=result.get_unitary(circ, decimals=3)
print(u1)
print(u2)
print(u1-u2)

pass_ = Unroller(['u1', 'u2', 'u3', 'cx'])
pm = PassManager(pass_)
new_circ = pm.run(circ)
print(new_circ)
