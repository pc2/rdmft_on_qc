import os
import copy
import math
import sys
import groundstate
import aca
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

try:
    from dmrgpy import fermionchain
except ImportError:
    print("please install https://github.com/joselado/dmrgpy")
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
                #print(op.primitive[i])
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
    return {"ops":ops,"qcs":qcs,"mesq":mesq,"coeff":coeff,"const":const}


dotenv.load_dotenv()
apikey=os.getenv("QISKIT_APIKEY")

tsim=True
L=8 #number of sites
Llocal=1 #number of sites in local approximation
shots=1024 #number of shots for qc-meaurements
seed=424242 #random seed






np.random.seed(seed)



#GS of Hubbard chain
t=-1.0
U=2.0
mu=-0.5*U
[E,D,W]=groundstate.hubbard_chain(L,t,U,mu,mode="DMRG")
print("E=",E)
print("D=",D)
print("W=",W)

#apply local approximation and ACA
print("local approximation with ",Llocal,"sites")
Wlocal=[]
W=[]
I=0
for i in range(L):
    if len(W)==Llocal:
        Wlocal.append(W)
        W=[]
    W.append(I)
    I=I+1
Wlocal.append(W)

print("local interactions on sites",Wlocal)

for i in range(len(Wlocal)):
    print("Functional for interaction on sites",Wlocal[i])
    orbinteract=[]
    
    #set up functional to be solved

    #reorder
    orbinteract=[]
    for j in Wlocal[i]:
        orbinteract.append(2*j)
        orbinteract.append(2*j+1)
    aca.printmat(2*L,2*L,"D_in_"+str(i),D)

    print("interacting orbitals:",orbinteract)

    Dreordered=aca.aca_reorder(2*L,orbinteract,D)
    aca.printmat(2*L,2*L,"D_reordered_"+str(i),Dreordered)

    Daca=Dreordered.copy()
    norb=2*L
    ninteract=len(orbinteract)
    for l in range(int(norb/ninteract-2)):
        Dacain=np.zeros((norb-l*ninteract,norb-l*ninteract),dtype=np.complex_)
        Dacain[::,::]=Daca[l*ninteract::,l*ninteract::]
        Dacaout=aca.aca(norb-l*ninteract,len(orbinteract),Dacain)
        Daca[l*ninteract::,l*ninteract::]=Dacaout
        aca.printmat(2*L,2*L,"D_aca_"+str(l),Daca)

#    Faca_local=hubbard_F(U,Wlocal,Daca)

    print("WARNING: This is completely UNoptimized code, for larger number of orbitals please use the Fortan implementation.")

    exit()

exit()





#qubit_converter = QubitConverter(mapper=ParityMapper())
quibit_converter = QubitConverter(mapper=JordanWignerMapper())
#qubit_converter = QubitConverter(mapper=BravyiKitaevMapper())

# setup the initial state for the ansatz
#entangler_map = [(0, 1), (1, 2), (2, 0)]
ansatz = TwoLocal(4,rotation_blocks = ['rx', 'ry'], entanglement_blocks = 'cz',entanglement='linear', reps=1, parameter_prefix = 'y')

print(ansatz)



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

#define registers
q = QuantumRegister(2*L)
c = ClassicalRegister(1)
qc=QuantumCircuit(q,c)
qc=qc.compose(ansatz)
qc=qc.decompose()


print("variational state:")
print(qc)

qc.draw(output='text',filename="ansatz.txt")
qc.draw(output='mpl',filename="ansatz.png")

qubits=[]
for q in range(2*L):
    qubits.append(q)

up=0
dn=1

c_ops=[]
for i in range(L):
    co=[]
    for j in range(2):
        co.append(FermionicOp("+_"+str(2*i+j),register_length=2*L))
    c_ops.append(co)
#build interaction operator 
interact=0
for i in range(L):
    interact+=U*(~c_ops[i][up] @ c_ops[i][up]@~c_ops[i][dn] @ c_ops[i][dn])
print(interact)

qc_interact=qcs_for_op(interact,qubit_converter,qc)
for i in range(len(qc_interact['qcs'])):
    print("op=",qc_interact['ops'][i],", measuring qubit",qc_interact['mesq'][i],", coeff=",qc_interact['coeff'][i])
    print(qc_interact['qcs'][i])

#build list of constraints
constraints=[]
for i in range(L):
    for si in range(2):
        for j in range(L):
            for sj in range(2):
                if i>j:
                    continue
                #real part
                c={}
                c['observable']='1RDM'
                c['type']='real'
                c['i']=i
                c['j']=j
                c['si']=si
                c['sj']=sj
                m=c_ops[i][si]@ ~ c_ops[j][sj]
                c['op']=0.5*(m+~m)
                c['qcs']=qcs_for_op(c['op'],qubit_converter,qc)
                c['cval']=0.5*(D[2*i+si,2*j+sj]+np.conjugate(D[2*i+si,2*j+sj])).real
                constraints.append(c)
                if i!=j:
                    c={}
                    c['observable']='1RDM'
                    c['type']='imag'
                    c['i']=i
                    c['j']=j
                    c['si']=si
                    c['sj']=sj
                    m=c_ops[i][si]@ ~ c_ops[j][sj]
                    c['op']=0.5/1j*(m-~m)
                    c['qcs']=qcs_for_op(c['op'],qubit_converter,qc)
                    c['cval']=(0.5/1j*(D[2*i+si,2*j+sj]-np.conjugate(D[2*i+si,2*j+sj]))).real
                    constraints.append(c)

print("number of constraints=",len(constraints))

c_qc=np.zeros(len(constraints))
c_exact=np.zeros(len(constraints))
W_qc=0
W_exact=0

qcs=[]
#programs for constraints

for i in range(len(constraints)):
    print(constraints[i])
    #evaluate constraint
    state=qc.bind_parameters(initial_point).decompose()
    if tsim:
        #check with exact value from expectation value of hermitian operator
        state.save_expectation_value(qubit_converter.convert(constraints[i]['op']),qubits)
        result=backend.run(state).result()
        exp = result.data()['expectation_value']
        c_exact[i]=exp
    for q in constraints[i]['qcs']['qcs']:
        q=q.bind_parameters(initial_point).decompose()
        qcs.append(q)

#programs for interaction        
if tsim:
    #check with exact value from expectation value of hermitian operator
    state=qc.bind_parameters(initial_point).decompose()
    state.save_expectation_value(qubit_converter.convert(interact),qubits)
    result=backend.run(state).result()
    exp = result.data()['expectation_value']
    W_exact=exp
for q in qc_interact['qcs']:
    q=q.bind_parameters(initial_point).decompose()
    qcs.append(q)
        

print("quantum programs=",len(qcs))        

#measure all constraints
jobs=execute(qcs,backend=backend,backend_properties=backend.properties(),shots=shots,seed_simulator=seed,seed_transpiler=seed)#,optimization_level=3)
if not tsim or True:
    print("waiting for job to finish: ",jobs.job_id())
    jobs.wait_for_final_state()
    print("job finished: ",jobs.job_id())


if jobs.done():
    res=jobs.result().results
    I=0
    #build constraint values from results
    for ic in range(len(constraints)):
        a=constraints[ic]['qcs']['const']
        for i in range(len(constraints[ic]['qcs']['qcs'])):
            v=-2*res[I].data.counts['0x1']/shots+1
            I=I+1
            a=a+v*constraints[ic]['qcs']['coeff'][i]
        c_qc[ic]=a
    #build interaction expectation value from result
    a=qc_interact['const']
    for i in range(len(qc_interact['qcs'])):
        v=-2*res[I].data.counts['0x1']/shots+1
        I=I+1
        a=a+v*qc_interact['coeff'][i]
    W_qc=a

print("W_exact=",W_exact)
print("W_qc=",W_qc)
for i in range(len(constraints)):
    print("c",constraints[i]['observable'],constraints[i]['type'],constraints[i]['i'],constraints[i]['si'],constraints[i]['j'],constraints[i]['sj'],"exact=",c_exact[i],"qc=",c_qc[i],"c=",constraints[i]['cval'])


#augmented Lagrangian
#for oiter in range(100):
    #unconstrainted steps
#    for iiter in range(100):
        #spsa step
    #multiplier and penalty update

