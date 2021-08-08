import os
import copy
import math
import sys
import groundstate
import aca
import ci
import graphclique
import time
from scipy import sparse
from scipy.optimize import minimize,BFGS
from math import pi
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
    from qiskit.circuit.library import TwoLocal,EfficientSU2
    #from qiskit.circuit import QuantumCircuit
    from qiskit.circuit import Parameter, ParameterVector, ParameterExpression
    from qiskit import execute
    from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
    from qiskit.algorithms import VQE
    from qiskit.algorithms.optimizers import L_BFGS_B,SPSA,COBYLA,QNSPSA
    from qiskit.opflow.primitive_ops import PauliOp
    from qiskit.opflow.state_fns import CircuitStateFn
    from qiskit.quantum_info import Pauli
    from qiskit.opflow.gradients import Gradient, NaturalGradient, QFI, Hessian
    from qiskit.opflow import Z, X, I, StateFn, CircuitStateFn, SummedOp
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

try:
    import configparser
except ImportError:
    print("installing configparser...")
    os.system("pip3 install configparser")
    print("installed configparser. Please restart.")
    quit()

def opt_callback(nfunc,par,f,stepsize,accepted):
    print("Opt step:",nfunc,par,f,stepsize,accepted)


def qcs_for_op(m,qubit_converter,qc,num_particles=0,tmeasure=True):
    m_op = qubit_converter.convert(m,num_particles=num_particles)
    ops=[]
    qcs=[]
    mesq=[]
    coeff=[]
    const=0

    #build measuring programs here
    p=m_op.to_pauli_op()
    #print(p)
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
            if tmeasure:
                q.measure(zs[0],0)
            qcs.append(q)
        else:
            const+=op.coeff
    return {"ops":ops,"qcs":qcs,"mesq":mesq,"coeff":coeff,"const":const}


dotenv.load_dotenv()
apikey=os.getenv("QISKIT_APIKEY")

config = configparser.ConfigParser()
config.sections()
config.read(sys.argv[1])

seed=int(config['rnd']['seed'])

two_qubit_reduction=config.getboolean('QC','two_qubit_reduction')
combine_qc_programs=config.getboolean('QC','combine_qc_programs')
tdoqc=config.getboolean('QC','tdoqc')
tsim=config.getboolean('QC','tsim')
tsampling=config.getboolean('QC','tsampling')
tnoise=config.getboolean('QC','tnoise')
if tnoise:
    tsampling=True

tcheck=config.getboolean('QC','tcheck')
tobj=config.getboolean('QC','tobj')
shots=int(config['QC']['shots'])
systemtype=config['system']['type']
L=int(config['system']['L'])
Lint=int(config['system']['Lint'])
if Lint==-1:
    Lint=L
Llocal=int(config['LocalApprox']['Llocal'])
ilocal=int(config['LocalApprox']['ilocal'])
np.random.seed(seed)

#GS of system to get a physical density matrix
if systemtype=='hubbard_chain':
    U=float(config['system']['U'])
    t=float(config['system']['t'])
    mu=float(config['system']['mu'])
    [E,D,W]=groundstate.hubbard_chain(L,t,U,mu,Lint,mode="DMRG")
else:
    print("system type not implemented")
    exit();

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

print("Functional for interaction on sites",Wlocal[ilocal])
orbinteract=[]

#set up functional to be solved
#reorder orbitals
norb=2*L

orbinteract=[]
ninteract=0
for j in Wlocal[ilocal]:
    ninteract+=2
    orbinteract.append(2*j)
    orbinteract.append(2*j+1)
aca.printmat(norb,norb,"D_in_"+str(ilocal),D)

print("interacting orbitals:",orbinteract)

[Dreordered,orbinteractaca]=aca.aca_reorder(norb,orbinteract,D,orbinteract)
aca.printmat(norb,norb,"D_reordered_"+str(ilocal),Dreordered)
print("orbinteractaca=",orbinteractaca)

Daca=Dreordered.copy()
Daca_levels=[]
#Daca_levels.append(np.copy(Daca[0:ninteract,0:ninteract]))

for l in range(0,int(norb/ninteract-1)):
    Dacain=np.zeros((norb-l*ninteract,norb-l*ninteract),dtype=np.complex_)
    Dacain[::,::]=np.copy(Daca[l*ninteract::,l*ninteract::])
    Dacaout=aca.aca(norb-l*ninteract,ninteract,Dacain)
    Daca[l*ninteract::,l*ninteract::]=np.copy(Dacaout)
    aca.printmat(norb,norb,"D_aca_"+str(l),Daca)
    aca.printmat((l+1)*ninteract,(l+1)*ninteract,"D_aca_truncated_"+str(l),Daca[0:(l+1)*ninteract,0:(l+1)*ninteract])
    Daca_levels.append(np.copy(Daca[0:(l+1)*ninteract,0:(l+1)*ninteract]))
Daca_levels.append(np.copy(Daca))

l=int(config['ACA']['level'])
Daca=Daca_levels[l]
norb_aca=np.shape(Daca)[0]
aca.printmat(norb_aca,norb_aca,"D_aca_used_",Daca)
print("aca: level",l," norb=",norb_aca)

print("Daca=",Daca)

#count spin-up and spin-down electrons
Naca1=0
Naca2=0
for i in range(int(norb_aca/2)):
    Naca1+=abs(Daca[2*i,2*i])
    Naca2+=abs(Daca[2*i+1,2*i+1])
if abs(Naca1-int(Naca1))<1e-4:
    Naca1=int(Naca1)
if abs(Naca2-int(Naca2))<1e-4:
    Naca2=int(Naca2)
Naca=(Naca1,Naca2)    
print("Naca=",Naca)

tdoci=config.getboolean('CI','tdoci')
if tdoci:
    tcplx=config.getboolean('CI','tcplx')
    options={'tol': float(config['CI']['tol']),'maxiter': int(config['CI']['maxiter'])}
    Faca_local=ci.F_hubbard(norb_aca,U,orbinteractaca,Daca,options,tcplx)
    print("F_aca (full space)=",Faca_local)

if tdoci:
    tcplx=config.getboolean('CI','tcplx')
    options={'tol': float(config['CI']['tol']),'maxiter': int(config['CI']['maxiter'])}
    Faca_local=ci.F_hubbard(norb_aca,U,orbinteractaca,Dreordered,options,tcplx)
    print("F_reordered (full space)=",Faca_local)

qubit_converter = QubitConverter(mapper=ParityMapper())
mapping=config['QC']['mapping']
if mapping=="Parity":
    qubit_converter = QubitConverter(mapper=ParityMapper(),two_qubit_reduction=two_qubit_reduction)
elif mapping=="JordanWigner":
    quibit_converter = QubitConverter(mapper=JordanWignerMapper(),two_qubit_reduction=two_qubit_reduction)
elif mapping=="BravyiKitaev":
    qubit_converter = QubitConverter(mapper=BravyiKitaevMapper(),two_qubit_reduction=two_qubit_reduction)
else:
    print("fermionic mapping unknown")
    exit()

if not tdoqc:
    exit()

entanglement=config['QC']['entanglement']
reps=int(config['QC']['reps'])
entanglement_blocks=config['QC']['entanglement_blocks']
rotation_blocks=config['QC']['rotation_blocks'].split(",")
if entanglement=="map":
    entanglement_map=config['QC']['entanglement_map']
    entanglement=[]
    for e in entanglement_map.split(","):
        entanglement.append((int(e.split("_")[0]),int(e.split("_")[1])))





# set the backend for the quantum computation
backend = BasicAer.get_backend('qasm_simulator')
optimization_level=1
if tnoise:
    backend = BasicAer.get_backend('qasm_simulator')
    optimization_level=1
else:
    #backend = BasicAer.get_backend('statevector_simulator')
    backend = Aer.get_backend('aer_simulator_statevector')
    optimization_level=0

backend_check = BasicAer.get_backend('statevector_simulator')
#backend_check = Aer.get_backend('aer_simulator_statevector')
if not tsim:
    IBMQ.save_account(apikey,overwrite=True)
    provider = IBMQ.load_account()
    print(provider.backends(n_qubits=5, operational=True))
    backend = provider.backend.ibmq_quito

num_particles=Naca

up=0
dn=1

c_ops=[]
for i in range(norb_aca):
    c_ops.append(FermionicOp("-_"+str(i),register_length=norb_aca))

#build interaction operator 
interact=0
Waca=[]
for i in range(int(len(orbinteract)/2)):
    Waca.append(i)

for i in Waca:
    if not tobj:
        #FIXME
        interact+=0.000001*((~c_ops[2*i]) @ c_ops[2*i]@(~c_ops[2*i+1]) @ c_ops[2*i+1])
    else:
        interact+=U*((~c_ops[2*i]) @ c_ops[2*i]@(~c_ops[2*i+1]) @ c_ops[2*i+1])


#get number of qubits
nq=norb_aca
if two_qubit_reduction:
    m_op = qubit_converter.convert(interact,num_particles=num_particles)
    nq=len(m_op.to_pauli_op()[0].primitive)
print("nq=",nq)    

# setup the initial state for the ansatz
ansatz = TwoLocal(nq,rotation_blocks = rotation_blocks, entanglement_blocks = entanglement_blocks,entanglement=entanglement, reps=reps, parameter_prefix = 'y',insert_barriers=False)
    
#define registers
q = QuantumRegister(nq)
c = ClassicalRegister(1)
qc=QuantumCircuit(q,c)
qc=qc.compose(ansatz)
qc=qc.decompose()

print("interact=",interact)
qc_interact=qcs_for_op(interact,qubit_converter,qc,num_particles)
for i in range(len(qc_interact['qcs'])):
    print("op=",qc_interact['ops'][i],", measuring qubit",qc_interact['mesq'][i],", coeff=",qc_interact['coeff'][i])
    print(qc_interact['qcs'][i])

iop=qubit_converter.convert(interact,num_particles=num_particles)
interact_sparse=sparse.csr_matrix(iop.to_matrix())

print("variational state:")
print(qc)

qc.draw(output='text',filename="ansatz.txt")
qc.draw(output='mpl',filename="ansatz.png")

if False:
    #gradient test
    q = QuantumRegister(nq)
    qc2=QuantumCircuit(q)
    qc2=qc2.compose(ansatz).decompose()
    #print(ansatz)
    H = qubit_converter.convert(interact,num_particles=num_particles) #qc_interact2['qcs'][i]
    params = ParameterVector('y', length=ansatz.num_parameters)
    print("params=",params)
    qc2=qc2.assign_parameters(params)

    op = ~StateFn(H) @ CircuitStateFn(primitive=qc2, coeff=1.)
    print("op=",op)


    grad = Gradient(grad_method='param_shift').convert(operator = op, params = params)
    for i in range(len(grad)):
        print("grad=",params[i],grad[i])
    exit()



qubits=[]
for q in range(nq):
    qubits.append(q)


#build list of constraints
constraints=[]
for i in range(norb_aca):
    for j in range(norb_aca):
        if i>j:
            continue
        #real part
        c={}
        c['observable']='1RDM'
        c['type']='real'
        c['i']=i
        c['j']=j
        m=~c_ops[i]@ c_ops[j]
        c['op']=0.5*(m+~m)
        op=qubit_converter.convert(c['op'],num_particles=num_particles)
        c['opsparse']=sparse.csr_matrix(op.to_matrix())
        c['qcs']=qcs_for_op(c['op'],qubit_converter,qc,num_particles)
        c['cval']=0.5*(Daca[i,j]+np.conjugate(Daca[i,j])).real
        constraints.append(c)
        if i!=j:
            c={}
            c['observable']='1RDM'
            c['type']='imag'
            c['i']=i
            c['j']=j
            m=~c_ops[i]@ c_ops[j]
            c['op']=0.5/1j*(m-~m)
            op=qubit_converter.convert(c['op'],num_particles=num_particles)
            c['opsparse']=sparse.csr_matrix(op.to_matrix())
            c['qcs']=qcs_for_op(c['op'],qubit_converter,qc,num_particles)
            c['cval']=(0.5/1j*(Daca[i,j]-np.conjugate(Daca[i,j]))).real
            constraints.append(c)

print("number of constraints=",len(constraints))

pauliops=[]
for cc in qc_interact['ops']:
    pauliops.append(str(cc.primitive))

for c in constraints:
    for cc in c['qcs']['ops']:
        pauliops.append(str(cc.primitive))

print("plain pauliops=",len(pauliops))
print(pauliops)       

pauliops2=[]
#get unique pauliops
for p in pauliops:
    found=False
    for p2 in pauliops2:
        if p==p2:
            found=True
            break
    if found:
        print("duplicate pauli op:",p)
    if not found:
        pauliops2.append(p)

print("unique pauliops=",len(pauliops2))
#print(pauliops2)       

if combine_qc_programs:
    #combine as many pauliops to quantum programs as possible
    nodes=[]
    for p in pauliops2:
        nodes.append(p)

    cliques=graphclique.cliquecover(nodes)
    exit()    

c_qc=np.zeros(len(constraints))
c_exact=np.zeros(len(constraints))
W_qc=0
W_exact=0

initial_point = np.random.random(ansatz.num_parameters)
print("number of parameters:",ansatz.num_parameters)

qcs=[]
#programs for constraints

circ=qc.bind_parameters(initial_point)
job=execute(circ,backend_check)
result=job.result()
psi=result.get_statevector()
for i in range(len(constraints)):
    print(constraints[i])
    #evaluate constraint
    if tcheck:
        exp=np.dot(np.conj(psi),constraints[i]['opsparse'].dot(psi)).real
        c_exact[i]=exp
    for q in constraints[i]['qcs']['qcs']:
        qcs.append(q)

#programs for interaction        
if tcheck:
    #check with exact value from expectation value of hermitian operator
    exp=np.dot(np.conj(psi),interact_sparse.dot(psi)).real
    #qc2=QuantumCircuit(QuantumRegister(nq)).compose(ansatz)
    #state=CircuitStateFn(qc2.bind_parameters(initial_point))
    #op=qubit_converter.convert(interact,num_particles=num_particles)
    #exp=state.adjoint().compose(op).compose(state).eval().real
    W_exact=exp
for q in qc_interact['qcs']:
    qcs.append(q)
        

print("quantum programs=",len(qcs))        


if tcheck:
#measure all constraints
    qcsp=[]
    for q in qcs:
        qcsp.append(q.bind_parameters(initial_point).decompose())
    jobs=execute(qcsp,backend=backend,backend_properties=backend.properties(),shots=shots,seed_simulator=seed,seed_transpiler=seed)#,optimization_level=3)
    if not tsim: 
        print("waiting for job to finish: ",jobs.job_id())
        jobs.wait_for_final_state(wait=1)
        print("job finished: ",jobs.job_id())
    else:
        jobs.wait_for_final_state(wait=0.05)

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
        print("c",constraints[i]['observable'],constraints[i]['type'],constraints[i]['i'],"exact=",c_exact[i],"qc=",c_qc[i],"c=",constraints[i]['cval'])

c_qc=np.zeros(len(constraints))

rdmf_obj_eval=0
rdmf_cons_eval=0
tprintevals=True


def rdmf_obj(x):
    global rdmf_obj_eval
    global tprintevals
    t1=time.perf_counter()
    t2=0
    t3=0
    t4=0
    t5=0
    rdmf_obj_eval+=1
    
    w_qc=0
    L=0
    if not tsampling:
        circ=qc.bind_parameters(x)
        job=execute(circ,backend_check)
        result=job.result()
        psi=result.get_statevector()
        t2=time.perf_counter()
        #build constraint values from results
        for ic in range(len(constraints)):
            a=np.dot(np.conj(psi),constraints[ic]['opsparse'].dot(psi)).real
            c_qc[ic]=a-constraints[ic]['cval']
            L=L+lagrange[ic]*c_qc[ic]+0.5*penalty*(c_qc[ic])**2
        t3=time.perf_counter()
        #build interaction expectation value from result
        a=np.dot(np.conj(psi),interact_sparse.dot(psi)).real
        t4=time.perf_counter()
        W_qc=a
        L=L+W_qc
        t5=time.perf_counter()
        
    else:
        qcsp=[]
        for q in qcs:
            qcsp.append(q.bind_parameters(x))#.decompose())
        t2=time.perf_counter()
        jobs=execute(qcsp,backend=backend,backend_properties=backend.properties(),shots=shots,seed_simulator=seed+rdmf_obj_eval,seed_transpiler=seed,optimization_level=optimization_level)
        t3=time.perf_counter()
        if not tsim: 
            print("waiting for job to finish: ",jobs.job_id())
            jobs.wait_for_final_state(wait=1)
            print("job finished: ",jobs.job_id())
        else:
            jobs.wait_for_final_state(wait=0.05)
        t4=time.perf_counter()

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
                c_qc[ic]=a-constraints[ic]['cval']
                L=L+lagrange[ic]*c_qc[ic]+0.5*penalty*(c_qc[ic])**2
            #build interaction expectation value from result
            a=qc_interact['const']
            for i in range(len(qc_interact['qcs'])):
                v=-2*res[I].data.counts['0x1']/shots+1
                I=I+1
                a=a+v*qc_interact['coeff'][i]
            W_qc=a
            L=L+W_qc
        t5=time.perf_counter()
    if tprintevals:
        print("L=",L,"W=",W_qc,"sum(c^2)=",np.sum(c_qc**2),"t=",t2-t1,t3-t2,t4-t3,t5-t4)
    return L

def rdmf_cons(x):
    global rdmf_cons_eval
    global tprintevals
    c_qc=np.zeros(len(constraints))
    t1=time.perf_counter()
    t2=0
    t3=0
    t4=0
    t5=0
    rdmf_cons_eval+=1
    
    w_qc=0
    L=0
    if not tsampling:
        circ=qc.bind_parameters(x)
        job=execute(circ,backend_check)
        result=job.result()
        psi=result.get_statevector()
        t2=time.perf_counter()
        #build constraint values from results
        for ic in range(len(constraints)):
            a=np.dot(np.conj(psi),constraints[ic]['opsparse'].dot(psi)).real
            c_qc[ic]=a-constraints[ic]['cval']
            L=L+lagrange[ic]*c_qc[ic]+0.5*penalty*(c_qc[ic])**2
        t3=time.perf_counter()
        #build interaction expectation value from result
        a=np.dot(np.conj(psi),interact_sparse.dot(psi)).real
        t4=time.perf_counter()
        W_qc=a
        L=L+W_qc
        t5=time.perf_counter()
        
    else:
        qcsp=[]
        for q in qcs:
            qcsp.append(q.bind_parameters(x))#.decompose())
        t2=time.perf_counter()
        jobs=execute(qcsp,backend=backend,backend_properties=backend.properties(),shots=shots,seed_simulator=seed+rdmf_obj_eval,seed_transpiler=seed,optimization_level=optimization_level)
        t3=time.perf_counter()
        if not tsim: 
            print("waiting for job to finish: ",jobs.job_id())
            jobs.wait_for_final_state(wait=1)
            print("job finished: ",jobs.job_id())
        else:
            jobs.wait_for_final_state(wait=0.05)
        t4=time.perf_counter()

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
                c_qc[ic]=a-constraints[ic]['cval']
                L=L+lagrange[ic]*c_qc[ic]+0.5*penalty*(c_qc[ic])**2
            #build interaction expectation value from result
            a=qc_interact['const']
            for i in range(len(qc_interact['qcs'])):
                v=-2*res[I].data.counts['0x1']/shots+1
                I=I+1
                a=a+v*qc_interact['coeff'][i]
            W_qc=a
            L=L+W_qc
        t5=time.perf_counter()
    if tprintevals:
        print("L=",L,"W=",W_qc,"sum(c^2)=",np.sum(c_qc**2),"t=",t2-t1,t3-t2,t4-t3,t5-t4)
    return c_qc

qiskit.utils.algorithm_globals.random_seed=seed
maxiter=int(config['QC']['maxiter'])    

x0=initial_point
tprintevals=False

if tsim and False:
    penalty=0
    lagrange=np.zeros(len(constraints))
    method="trust-constr"
#    method="SLSQP"
    print("minimizing with constrained minimization with "+method+" (no sampling, no noise)")
    eq_cons={'type': 'eq','fun' : lambda x: rdmf_cons(x)}
    if method=="trust-constr":
        res=minimize(rdmf_obj, x0, method=method, constraints=[eq_cons],tol=1e-4,options={'maxiter':10000,'verbose': 3,'disp': True,'initial_tr_radius':4*pi,'gtol':1e-4,'xtol':1e-4})
    else:
        res=minimize(rdmf_obj, x0, method=method, constraints=[eq_cons],tol=1e-4,options={'maxiter':10000,'verbose': 3,'iprint':2,'disp': True})
    print(res)

penalty=10
lagrange=np.zeros(len(constraints))
#augmented Lagrangian
for oiter in range(100):
    value=0
    point=[]
    nfev=0
    #unconstrainted steps
    if tsim and not tsampling and not tnoise:
        print("minimizing with augmented Lagrangian and LBFGS_B (no sampling, no noise)")
        method="BFGS"
        tprintevals=True
        res=minimize(rdmf_obj, x0, method=method,tol=1e-2,options={'maxiter':10000,'verbose': 2,'disp': True})
        point=res.x
        value=rdmf_obj(point)
        nfev=res.nfev
        #optimizer = L_BFGS_B(maxiter=maxiter,tol=1e-3)
        #[point, value, nfev]=optimizer.optimize(num_vars=ansatz.num_parameters,objective_function=rdmf_obj,initial_point=x0,approx_grad=True)
    elif tsim and tsampling and not tnoise:
        print("minimizing with augmented Lagrangian and COBYLA (with sampling, no noise)")
        tprintevals=True
        optimizer = COBYLA(maxiter=maxiter,disp=True,tol=1e-2)
        [point, value, nfev]=optimizer.optimize(num_vars=ansatz.num_parameters,objective_function=rdmf_obj,initial_point=x0)
    elif not tsim:
        print("minimizing with augmented Lagrangian and SPSA (with sampling, with noise)")
        tprintevals=True
        optimizer = SPSA(maxiter=maxiter,second_order=False)#,callback=opt_callback)
        print("stddev(L)=",optimizer.estimate_stddev(rdmf_obj,initial_point,avg=25))
        print("calibrating")
        [learning_rate,perturbation]=optimizer.calibrate(rdmf_obj,initial_point,stability_constant=0, target_magnitude=None, alpha=0.602, gamma=0.101, modelspace=False)
        print("learning_rate=",learning_rate)
        print("perturbation=",perturbation)
        optimizer = SPSA(maxiter=maxiter,second_order=False,perturbation=perturbation,learning_rate=learning_rate)#,callback=opt_callback)
        print("minimizing")
        [point, value, nfev]=optimizer.optimize(num_vars=ansatz.num_parameters,objective_function=rdmf_obj,initial_point=initial_point)

    print("point=",point)
    print("value=",value)
    print("nfev=",nfev)
    print("constraint violation sum(c^2)=",np.sum(c_qc**2))

    x0=point
    #multiplier and penalty update
    for i in range(len(constraints)):
        lagrange[i]=lagrange[i]+penalty*c_qc[i]
    penalty=penalty*2
    print("lagrange=",lagrange)
    print("penalty=",penalty)




