import numpy as np
from qiskit import *
from qiskit import Aer
from qiskit.circuit import Gate
from math import pi
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Unroller
from qiskit.quantum_info import Pauli
from qiskit.opflow.primitive_ops import PauliOp

def cnot(nq,c,t):
#Changing the simulator 
    backend = Aer.get_backend('unitary_simulator')

#The circuit without measurement
    circ = QuantumCircuit(nq)
    circ.cx(c,t)
    job = execute(circ, backend)
    result = job.result()

    pass_ = Unroller(['u1', 'u2', 'u3', 'cx'])
    pm = PassManager(pass_)
    new_circ = pm.run(circ)
    return result.get_unitary(circ, decimals=3)

def h(nq,c):
#Changing the simulator 
    backend = Aer.get_backend('unitary_simulator')

#The circuit without measurement
    circ = QuantumCircuit(nq)
    circ.h(c)
    job = execute(circ, backend)
    result = job.result()

    pass_ = Unroller(['u1', 'u2', 'u3', 'cx'])
    pm = PassManager(pass_)
    new_circ = pm.run(circ)
    return result.get_unitary(circ, decimals=3)


nq=4
iiiz=Pauli('IIIZ')
izzz=Pauli('IZZZ')
zizi=Pauli('ZIZI')
izii=Pauli('IZII')
opin=[iiiz,izzz,zizi,izii]

ziii=Pauli('ZIII')
iizi=Pauli('IIZI')
op2=[ziii,izii,iizi,iiiz]

U=cnot(nq,1,3)
U=np.matmul(cnot(nq,0,1),U)
U=np.matmul(cnot(nq,2,1),U)

for op in opin:
    print(op)
    t=np.matmul(np.matmul(U,op.to_matrix()),np.transpose(U))
    for opref in op2:
        print(np.sum(abs(t-opref.to_matrix())))


nq=1
x=Pauli('X')
y=Pauli('Y')
z=Pauli('Z')
U=h(nq,0)

for op in [x,y,z]:
    print(op)
    t=np.matmul(np.matmul(U,op.to_matrix()),np.transpose(U))
    for opref in [x,y,z]:
        print(np.sum(abs(t-opref.to_matrix())))



