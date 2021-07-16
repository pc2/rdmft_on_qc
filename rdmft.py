import os
import copy
try:
    import numpy as np
except ImportError:
    print("installing numpy...")
    os.system("pip3 install numpy")
    print("installed numpy. Please restart.")
    quit()


try:
    import dotenv
except ImportError:
    print("installing python-dotenv...")
    os.system("pip3 install python-dotenv")
    print("installed python-dotenv. Please restart.")
    quit()


try:
    import matplotlib.pyplot as plt
except ImportError:
    print("installing matplotlib...")
    os.system("pip3 install matplotlib")
    print("installed matplotlib. Please restart.")
    quit()

try:
    import qiskit
    from qiskit.circuit.library import TwoLocal
    #from qiskit.circuit import QuantumCircuit
    from qiskit import execute
    from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
    from qiskit.algorithms import VQE
    from qiskit.algorithms.optimizers import L_BFGS_B,SPSA
    from qiskit.opflow.primitive_ops import PauliOp
    from qiskit.quantum_info import Pauli
except ImportError:
    print("installing qiskit...")
    os.system("pip3 install qiskit")
    print("installed qiskit. Please restart.")
    quit()

try:
    import pyscf
except ImportError:
    print("installing pyscf...")
    os.system("pip3 install pyscf")
    print("installed pyscf. Please restart.")
    quit()

try:
    import qiskit_nature
    from qiskit import IBMQ, assemble, transpile,Aer
    from qiskit import BasicAer
    from qiskit_nature.algorithms import GroundStateEigensolver
    from qiskit_nature.converters.second_quantization import QubitConverter
    from qiskit_nature.drivers import PySCFDriver, UnitsType, Molecule
    from qiskit_nature.mappers.second_quantization import ParityMapper
    from qiskit_nature.problems.second_quantization import ElectronicStructureProblem
    from qiskit_nature.problems.second_quantization.electronic import ElectronicStructureProblem
    from qiskit_nature.converters.second_quantization import QubitConverter
    from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper, BravyiKitaevMapper    
    from qiskit_nature.operators.second_quantization import FermionicOp
    
except ImportError  as error:
    print("installing qiskit-nature because",error.message)
    os.system("pip3 install qiskit-nature")
    print("installed qiskit-nature. Please restart.")
    quit()

try:
    import pylatexenc
except ImportError  as error:
    print("installing pylatexenc")
    os.system("pip3 install pylatexenc")
    print("installed pylatexenc. Please restart.")
    quit()

try:
    import sympy
except ImportError:
    print("installing sympy...")
    os.system("pip3 install sympy")
    print("installed sympy. Please restart.")
    quit()

def opt_callback(nfunc,par,f,stepsize,accepted):
    print("Opt step:",nfunc,par,f,stepsize,accepted)

def qcs_for_op(m,qubit_converter,qc):
    m_op = qubit_converter.convert(m) #,num_particles=num_particles)
    ops=[]
    qcs=[]
    mesq=[]
    coeff=[]
    const=0

    #build measuring programs here
    p=m_op.to_pauli_op()
    for op in p:
        if op.to_circuit().depth()>0:
            coeff.append(op.coeff)
            ops.append(op)
            #map to only sigma_z
            q=copy.deepcopy(qc)

            sg=[]
            for i in range(len(op.primitive)):
                print(op.primitive[i])
                if str(op.primitive[i])=='Z':
                    sg.append("Z")            
                elif str(op.primitive[i])=='X':
                    q.h(i)
                    sg.append("Z")            
                elif str(op.primitive[i])=='Y':
                    q.sdg(i)
                    q.h(i)
                    sg.append("Z")            
                else:
                    sg.append("I")            
            #map multiple sigma_z to single sigma_z
            zs=[i for i,x in enumerate(sg) if x=='Z']
            mesq.append(zs[0])
            for i in range(1,len(zs)):
                q.cx(zs[i],zs[0])
            q.measure(zs[0],0)
            qcs.append(q)
        else:
            const+=op.coeff
    return [ops,qcs,mesq,coeff,const]



dotenv.load_dotenv()
apikey=os.getenv("QISKIT_APIKEY")

tsim=False
L=2
shots=1024
seed=424242

#VQE in qiskit from https://github.com/Qiskit/qiskit-nature/
np.random.seed(seed)

# Use PySCF, a classical computational chemistry software
# package, to compute the one-body and two-body integrals in
# electronic-orbital basis, necessary to form the Fermionic operator

#molecule = Molecule(geometry=[['H', [0., 0., 0.]],['H', [0., 0., 0.735]]],charge=0, multiplicity=1)
#driver = PySCFDriver(molecule=molecule,unit=UnitsType.ANGSTROM,basis='sto3g')
#problem = ElectronicStructureProblem(driver)

# generate the second-quantized operators
#second_q_ops = problem.second_q_ops()
#print("second_q_ops=",second_q_ops)

#main_op = second_q_ops[0]

#num_particles = (problem.molecule_data_transformed.num_alpha,problem.molecule_data_transformed.num_beta)
#num_spin_orbitals = 2 * problem.molecule_data.num_molecular_orbitals

#qubit_converter = QubitConverter(mapper=ParityMapper())
#qubit_converter = QubitConverter(mapper=JordanWignerMapper())
qubit_converter = QubitConverter(mapper=BravyiKitaevMapper())

# setup the initial state for the ansatz
#entangler_map = [(0, 1), (1, 2), (2, 0)]
ansatz = TwoLocal(4,rotation_blocks = ['rx', 'ry'], entanglement_blocks = 'cz',entanglement='linear', reps=1, parameter_prefix = 'y')

print(ansatz)

ansatz.draw(output='text',filename="ansatz.txt")
ansatz.draw(output='mpl',filename="ansatz.png")


# set the backend for the quantum computation
backend = Aer.get_backend('aer_simulator_statevector')
#backend = BasicAer.get_backend('qasm_simulator')
if not tsim:
    IBMQ.save_account(apikey,overwrite=True)
    provider = IBMQ.load_account()
    print(provider.backends(n_qubits=5, operational=True))
    backend = provider.backend.ibmq_quito

#optimizer = L_BFGS_B()
#optimizer = SPSA(maxiter=100,callback=opt_callback)
optimizer = SPSA(maxiter=100)
print(optimizer.print_options())
initial_point = np.random.random(ansatz.num_parameters)

if tsim and False:
    # setup and run VQE
    algorithm = VQE(ansatz=ansatz,optimizer=optimizer,quantum_instance=backend,initial_point=initial_point)
    calc = GroundStateEigensolver(qubit_converter, algorithm)
    res = calc.solve(problem)
    print(res)

#
q = QuantumRegister(2*L)
c = ClassicalRegister(1)
qc=QuantumCircuit(q,c)
qc=qc.compose(ansatz)

qc=qc.bind_parameters(initial_point).decompose()
print(qc)

qubits=[]
for q in range(2*L):
    qubits.append(q)

m=FermionicOp("+_0",register_length=2*L) @ FermionicOp("-_3",register_length=2*L)
m=m+~m
#m=FermionicOp("+_0",register_length=2*L) @ FermionicOp("-_0",register_length=2*L)@FermionicOp("+_1",register_length=2*L) @ FermionicOp("-_1",register_length=2*L)

[ops,qcs,mesq,coeff,const]=qcs_for_op(m,qubit_converter,qc)


for i in range(len(qcs)):
    print("op=",ops[i],", measuring qubit",mesq[i],", coeff=",coeff[i])
    print(qcs[i])

#measure everything
#check with exact value from expectation value of hermitian operator
if tsim:
    qc.save_expectation_value(qubit_converter.convert(m),[0,1,2,3])
    result=backend.run(qc).result()
    exp = result.data()['expectation_value']
    print("<psi|Op|psi>=",exp)


#tqcs=transpile(qcs,backend,seed_transpiler=seed,optimization_level=3)
#qobj = assemble(tqcs,backend)

jobs=execute(qcs,backend=backend,backend_properties=backend.properties(),shots=shots,seed_simulator=seed,seed_transpiler=seed)#,optimization_level=3)
print("waiting for job to finish: ",jobs.job_id())
jobs.wait_for_final_state()
print("job finished: ",jobs.job_id())

if jobs.done():
    res=jobs.result().results
    print(res)
    a=const
    for i in range(len(res)):
        print(ops[i]," --> ",res[i].data.counts)
        v=-2*res[i].data.counts['0x1']/shots+1
        a=a+v*coeff[i]
    print("a=",a)



exit()


#build Hubbard model
c_ops=[]
L=4
for i in range(L):
    c_site=[]
    for s in range(2):
        c_site.append(FermionicOp("-_"+str(2*i+s),register_length=2*L))
    c_ops.append(c_site)

print(c_ops)

hopping=0 #*FermionicOp("I_0",register_length=2*L)
for i in range(L-1):
    for s in range(2):
        hopping+=~c_ops[i][s] @ c_ops[i+1][s]
hopping+=~hopping
print(hopping)

interact=0
for i in range(L):
    interact+=~c_ops[i][0] @ c_ops[i][0]@~c_ops[i][1] @ c_ops[i][1]
print(interact)

#build 1-RDM qubit op
RDM=0

#build ee-interaction qubit op


#read in density matrix

#augmented Lagrangian
