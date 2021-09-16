import os
import copy
import math
import sys
import groundstate
import aca
import ci
import time
import measurement_circuits
import numpy as np
import dotenv
import matplotlib.pyplot as plt
import pylatexenc
import configparser
import itertools
import qiskit
from scipy import sparse
from scipy.optimize import minimize, BFGS
from math import pi
from qiskit.circuit.library import TwoLocal, EfficientSU2
from qiskit.circuit import Parameter, ParameterVector, ParameterExpression
from qiskit import IBMQ, assemble, transpile, Aer, BasicAer, execute, QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import L_BFGS_B, SPSA, COBYLA, QNSPSA
from qiskit.opflow.primitive_ops import PauliOp
from qiskit.opflow.state_fns import CircuitStateFn
from qiskit.quantum_info import Pauli
from qiskit.opflow.gradients import Gradient, NaturalGradient, QFI, Hessian
from qiskit.opflow import Z, X, I, StateFn, CircuitStateFn, SummedOp
from qiskit.providers.aer.noise import NoiseModel
from qiskit.visualization import plot_histogram
import qiskit_nature
from qiskit_nature.circuit.library import HartreeFock
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import ParityMapper
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper, BravyiKitaevMapper
from qiskit_nature.operators.second_quantization import FermionicOp

def opt_callback(nfunc,par,f,stepsize,accepted):
    #print("Opt step:",nfunc,par,f,stepsize,accepted)
    print("Opt step:",nfunc,f,stepsize,accepted)

def paulis_for_op(m,qubit_converter,qc,num_particles=0,tmeasure=True):
    m_op = qubit_converter.convert(m,num_particles=num_particles)
    pauliops=[]
    coeff=[]
    const=0

    #build measuring programs here
    p=m_op.to_pauli_op()
    for op in p:
        if op.to_circuit().depth()>0:
            coeff.append(op.coeff)
            pauliops.append(op)
        else:
            const+=op.coeff
    return {"pauliops":pauliops,"coeff":coeff,"const":const}


def measure_all_programs(nq,ansatz,mqcs,x,backend,shots,seed,tsim,optimization_level=3):
    #measure all quantum programs
    qcs_par=[]
    for mqc in mqcs:
        q = QuantumRegister(nq)
        c = ClassicalRegister(nq)
        qc=QuantumCircuit(q,c)
        qc=qc.compose(ansatz)
        qc=qc.compose(mqc['mqc'])
        qcs_par.append(qc.bind_parameters(x))
    if not tnoise:
        jobs = execute(
                qcs_par,
                backend=backend,
                backend_properties=backend.properties(), # why not use backend_ibmq_manila.properties() here?
                shots=shots,
                seed_simulator=seed,
                seed_transpiler=seed,
                optimization_level=optimization_level
            )
    else:
        #print("EXECUTING WITH NOISE!")
        if tsim:
            jobs = execute(
                qcs_par,
                backend=backend,
                backend_properties=ibmq_backend_props,  # ibmq_manila or simulator backend?
                shots=shots, # ?
                seed_simulator=seed, # ?
                seed_transpiler=seed, # ?
                optimization_level=optimization_level, # ?,
                coupling_map=ibmq_coupling_map,
                basis_gates=ibmq_basis_gates,
                noise_model=ibmq_noise_model
            )
        else:
            jobs = execute(
                qcs_par,
                backend=backend,
                backend_properties=ibmq_backend_props,  # ibmq_manila or simulator backend?
                shots=shots, # ?
                seed_simulator=seed, # ?
                seed_transpiler=seed, # ?
                optimization_level=optimization_level
            )
    
    if not tsim: 
        print("waiting for job to finish: ",jobs.job_id())
        jobs.wait_for_final_state(wait=1)
        print("job finished: ",jobs.job_id())
    else:
        jobs.wait_for_final_state(wait=0.01)

    if jobs.done():
        res=jobs.result().results
        for j in range(len(mqcs)):
            counts = jobs.result().get_counts(j)
            #print("NOISE count %i:" % j)
            #plot_histogram(counts).savefig("counts_%i.png" % j)

        Vs=[]
        #loop over programs
        for im in range(len(mqcs)):
            #for r in res[im].data.counts:
            #    print(r,bin(int(r, base=16))[2:].zfill(nq),res[im].data.counts[r]/shots)

            Va=[]
            #loop over measurements in program
            for ic in range(len(mqcs[im]['ops'])):
                #print("op=",mqcs[im]['ops'][ic],"is measured at",mqcs[im]['mqubits'][ic],mqcs[im]['signs'][ic])
                V=0
                for r in res[im].data.counts:
                    b=bin(int(r, base=16))[2:].zfill(nq)
                    v=res[im].data.counts[r]/shots
                    #if b[nq-1-mqcs[im]['mqubits'][ic]]=='1':
                    if b[nq-1-ic]=='1':
                        V=V-v
                    else:
                        V=V+v
                V=V*mqcs[im]['signs'][ic]
                Va.append(V)
            Vs.append(Va)
            #print(mqcs[im]['ops'],Va)
    return Vs

def measurements_to_interact(mqcs,v,pauli_interact):
    W=pauli_interact['const']

    for ic in range(len(pauli_interact['pauliops'])):
        op=str(pauli_interact['pauliops'][ic].primitive)
        #find measurement in mqcs
        found=False
        for im in range(len(mqcs)):
            if found:
                break
            for io in range(len(mqcs[im]['ops'])):
                if op==mqcs[im]['ops'][io]:
                    W=W+v[im][io]*pauli_interact['coeff'][ic]
                    found=True
                    break
    return W

def measurements_to_constraints(mqcs,v,constraints):
    c=np.zeros(len(constraints))
    for ic in range(len(constraints)):
        cv=constraints[ic]['pauli']['const']-constraints[ic]['cval']
        for ip in range(len(constraints[ic]['pauli']['pauliops'])):
            found=False
            op=str(constraints[ic]['pauli']['pauliops'][ip].primitive)
            #find measurement in mqcs
            for im in range(len(mqcs)):
                if found:
                    break
                for io in range(len(mqcs[im]['ops'])):
                    if op==mqcs[im]['ops'][io]:
                        cv=cv+v[im][io]*constraints[ic]['pauli']['coeff'][ip]
                        found=True
                        break
        c[ic]=cv
    return c

def rdmf_obj(x):
    global rdmf_obj_eval
    rdmf_obj_eval+=1
    c_qc=np.zeros(len(constraints))
    W_qc=0
    L=0

    if texact_expect:
#    if True:
        circ=qc.bind_parameters(x)
        job=execute(circ,backend_check)
        result=job.result()
        psi=result.get_statevector()
        #build constraint values from results
        for ic in range(len(constraints)):
            if not build_sparse:
                raise RuntimeError("build_sparse not enabled")
            a=np.dot(np.conj(psi),constraints[ic]['opsparse'].dot(psi)).real
            c_qc[ic]=a-constraints[ic]['cval']
            L=L+lagrange[ic]*c_qc[ic]+0.5*penalty*abs(c_qc[ic])**penalty_exponent
        #build interaction expectation value from result
        if not build_sparse:
            raise RuntimeError("build_sparse not enabled")
        a=np.dot(np.conj(psi),interact_sparse.dot(psi)).real
        W_qc=a
        L=L+W_qc
        L2=L
    else:
        L=0
        v=measure_all_programs(nq,ansatz,mqcs,x,backend,shots,seed+rdmf_obj_eval,tsim,optimization_level)
        W_qc=measurements_to_interact(mqcs,v,pauli_interact)
        c_qc=measurements_to_constraints(mqcs,v,constraints)
        for ic in range(len(constraints)):
            L=L+lagrange[ic]*c_qc[ic]+0.5*penalty*abs(c_qc[ic])**penalty_exponent
        L=L+W_qc
    if tprintevals and rdmf_obj_eval%printevalsevery==0:
        print("L=",L,"W=",W_qc,"sum(c^2)=",np.sum(c_qc**2))
    return L

def rdmf_cons(x):
    global rdmf_cons_eval
    c_qc=np.zeros(len(constraints))
    rdmf_cons_eval+=1
    
    L=0
    if texact_expect:
        circ=qc.bind_parameters(x)
        job=execute(circ,backend_check)
        result=job.result()
        psi=result.get_statevector()
        #build constraint values from results
        for ic in range(len(constraints)):
            if not build_sparse:
                raise RuntimeError("build_sparse not enabled")
            a=np.dot(np.conj(psi),constraints[ic]['opsparse'].dot(psi)).real
            c_qc[ic]=a-constraints[ic]['cval']
    else:
        raise RuntimeError("rdmft_cons must not be used with tsampling=True or tsim=False")
    return c_qc

#main code

dotenv.load_dotenv()
apikey=os.getenv("QISKIT_APIKEY")

config = configparser.ConfigParser()
config.sections()
config.read(sys.argv[1])

seed=int(config['rnd']['seed'])

build_sparse=config.getboolean('QC','build_sparse')
commutation_mode=config['QC']['commutation']
two_qubit_reduction=config.getboolean('QC','two_qubit_reduction')
tdoqc=config.getboolean('QC','tdoqc')
tsim=config.getboolean('QC','tsim')
tsimcons=config.getboolean('QC','tsimcons')
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
elif systemtype=='hubbard_2d':
    U=float(config['system']['U'])
    t=float(config['system']['t'])
    mu=float(config['system']['mu'])
    [E,D,W]=groundstate.hubbard_2d(L,t,U,mu,Lint,mode="DMRG")
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
if l==-1:
    l=len(Daca_levels)-1
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

if abs(Naca1-round(Naca1))<1e-4:
    print("rounding Naca1 to integer")
    Naca1=round(Naca1)
if abs(Naca2-round(Naca2))<1e-4:
    print("rounding Naca2 to integer")
    Naca2=round(Naca2)
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
    qubit_converter = QubitConverter(mapper=JordanWignerMapper(),two_qubit_reduction=two_qubit_reduction)
elif mapping=="BravyiKitaev":
    qubit_converter = QubitConverter(mapper=BravyiKitaevMapper(),two_qubit_reduction=two_qubit_reduction)
else:
    print("fermionic mapping unknown")
    exit()




# set the backend for the quantum computation
if tnoise:
    # Build noise model from backend properties
    provider = IBMQ.load_account()
    # >>> [b.name() for b in provider.backends()]
    # ['ibmq_qasm_simulator', 'ibmq_armonk', 'ibmq_santiago', 'ibmq_bogota', 'ibmq_lima', 'ibmq_belem', 'ibmq_quito', 'simulator_statevector', 'simulator_mps', 'simulator_extended_stabilizer', 'simulator_stabilizer', 'ibmq_manila']

    ibmq_backend = provider.get_backend(config['QC']['qc'])
    ibmq_backend_config = ibmq_backend.configuration()

    if ibmq_backend_config.n_qubits < norb_aca:
        raise RuntimeError("chosen qc doesn't have enough qubits")

    # >>> backend_config.basis_gates
    # ['id', 'rz', 'sx', 'x', 'cx', 'reset']
    ibmq_coupling_map = ibmq_backend_config.coupling_map

    #set coupling map for transpiler that is used in the construction of measurement circuits
    couplingstring=""
    for ic in range(len(ibmq_coupling_map)-1):
        couplingstring+=str(ibmq_coupling_map[ic][0])+"_"+str(ibmq_coupling_map[ic][1])+","
    ic=len(ibmq_coupling_map)-1
    couplingstring+=str(ibmq_coupling_map[ic][0])+"_"+str(ibmq_coupling_map[ic][1])
    config['QC']['transpiler_couplings']=couplingstring
    print("setting transpiler_couplings to",couplingstring)
    
    #set gates for transpiler that is used in the construction of measurement circuits
    config['QC']['transpiler_gates']=",".join(ibmq_backend_config.basis_gates)
    print("setting transpiler_gates to",",".join(ibmq_backend_config.basis_gates))

    # [[0, 1], [1, 0], [1, 2], [1, 3], [2, 1], [3, 1], [3, 4], [4, 3]]
    # Below: directed graph specifying fixed coupling. Nodes correspond to physical qubits (integers) and directed edges correspond to permitted CNOT gates.
    # cm = qiskit.transpiler.CouplingMap(coupling_map)
    # >>> cm.distance_matrix
    # array([[0., 1., 2., 2., 3.],
    #        [1., 0., 1., 1., 2.],
    #        [2., 1., 0., 2., 3.],
    #        [2., 1., 2., 0., 1.],
    #        [3., 2., 3., 1., 0.]])
    ibmq_backend_props = ibmq_backend.properties()
    # nqubits = len(backend_props.qubits)

    # Noise....
    ibmq_noise_model = NoiseModel.from_backend(ibmq_backend)
    ibmq_basis_gates = ibmq_noise_model.basis_gates
    # ibmq_quito: ['cx', 'id', 'reset', 'rz', 'sx', 'x'] ('cx' ist the only two-qubit gate)

    #backend = BasicAer.get_backend('qasm_simulator')
    backend = Aer.get_backend('qasm_simulator')
    optimization_level=2
else:
    #backend = BasicAer.get_backend('statevector_simulator')
    backend = Aer.get_backend('aer_simulator_statevector')
    optimization_level=2

#set entangler map in hardware-effcieint trial state
entanglement=config['QC']['entanglement']
reps=int(config['QC']['reps'])
entanglement_blocks=config['QC']['entanglement_blocks']
rotation_blocks=config['QC']['rotation_blocks'].split(",")
if entanglement=="map":
    entanglement_map=config['QC']['entanglement_map']
    entanglement=[]
    for e in entanglement_map.split(","):
        entanglement.append((int(e.split("_")[0]),int(e.split("_")[1])))

if entanglement=="qc":
    #get entangler from hardware coupling map
    if not tnoise:
        raise RuntimeError("entanglement=qc requires tnoise=True")
    entanglement=[]
    for c in ibmq_coupling_map:
        if [c[1],c[0]] in entanglement:
            continue
        if c[0]<norb_aca and c[1]<norb_aca:
            entanglement.append([c[0],c[1]])
    print("using entangler=",entanglement)


backend_check = BasicAer.get_backend('statevector_simulator')
#backend_check = Aer.get_backend('aer_simulator_statevector')

if not tsim:
    #IBMQ.save_account(apikey,overwrite=True)
    provider = IBMQ.load_account()
    backend = provider.get_backend(config['QC']['qc'])
    backend_config = ibmq_backend.configuration()
    status = backend.status()
    
    if backend_config.n_qubits < norb_aca:
        raise RuntimeError("chosen qc doesn't have enough qubits")
    if not status.operational:
        raise RuntimeError("chosen qc is not operational")

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
#ansatz = HartreeFock(nq,Naca,qubit_converter)

#define registers
q = QuantumRegister(nq)
c = ClassicalRegister(nq)
qc=QuantumCircuit(q,c)
qc=qc.compose(ansatz)
qc=qc.decompose()

print("variational state:")
print(qc)

print("ansatz-state is written to ansatz.txt")
qc.draw(output='text',filename="ansatz.txt")
#qc.draw(output='mpl',filename="ansatz.png")

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

pauli_interact=paulis_for_op(interact,qubit_converter,qc,num_particles)
interact_op=qubit_converter.convert(interact,num_particles=num_particles)
if build_sparse:
    interact_sparse=sparse.csr_matrix(interact_op.to_matrix())

#build list of constraints
print("building constraints...")
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
        if build_sparse:
            c['opsparse']=sparse.csr_matrix(op.to_matrix())
        c['pauli']=paulis_for_op(c['op'],qubit_converter,qc,num_particles)
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
            if build_sparse:
                c['opsparse']=sparse.csr_matrix(op.to_matrix())
            c['pauli']=paulis_for_op(c['op'],qubit_converter,qc,num_particles)
            c['cval']=(0.5/1j*(Daca[i,j]-np.conjugate(Daca[i,j]))).real
            constraints.append(c)
        #print("constraint",i,j,c['pauli'])

print("number of constraints=",len(constraints))

############################################################################
#build measurement programs for 1rdm
pauliops=[]
for c in constraints:
    for cc in c['pauli']['pauliops']:
        pauliops.append(str(cc.primitive))
unique_pauliops=measurement_circuits.get_unique_pauliops(pauliops)['unique pauliops']
print("unique pauliops",len(unique_pauliops))
#compute commuting cliques
cliques=measurement_circuits.build_cliques(commutation_mode,unique_pauliops)
if config.getboolean('QC','end_after_cliques'):
    exit()
#build measurement circuits for ciques
print("constructing measurement programs for 1rdm with commutativity mode",commutation_mode)
mqcs_1rdm=measurement_circuits.build_measurement_circuits(cliques,commutation_mode,nq,config)

############################################################################
#build measurement programs for interaction
pauliops=[]
for cc in pauli_interact['pauliops']:
    pauliops.append(str(cc.primitive))
unique_pauliops=measurement_circuits.get_unique_pauliops(pauliops)['unique pauliops']
#compute commuting cliques
cliques=measurement_circuits.build_cliques(commutation_mode,unique_pauliops)
#build measurement circuits for ciques
print("constructing measurement programs for interaction with commutativity mode",commutation_mode)
mqcs_interact=measurement_circuits.build_measurement_circuits(cliques,commutation_mode,nq,config)

#############################################################################
#combine measurement programs for 1rdm and interaction
mqcs=copy.deepcopy(mqcs_1rdm)
mqcs.extend(mqcs_interact)

print("Number of quantum programs:",len(mqcs_1rdm),"for 1rdm and ",len(mqcs_interact),"for interaction")
print(mqcs)


#set intial parameters
np.random.seed(seed)
initial_point = np.random.random(ansatz.num_parameters)
print("number of parameters:",ansatz.num_parameters)

if not tdoqc:
    exit()

v=measure_all_programs(nq,ansatz,mqcs,initial_point,backend,shots,seed,tsim,optimization_level)
W_qc=measurements_to_interact(mqcs,v,pauli_interact)
c_qc=measurements_to_constraints(mqcs,v,constraints)

if tcheck:
    circ=ansatz.bind_parameters(initial_point)
    job=execute(circ,backend_check)
    result=job.result()
    counts = result.get_counts(0)
    #plot_histogram(counts).savefig("tcheck_counts.png")

    psi=result.get_statevector()
    #print("NONOISE STATE VECTOR:", psi)

    W_exact=np.dot(np.conj(psi),interact_sparse.dot(psi)).real
    print("W_exact=",W_exact,"W_qc=",W_qc)

    c_exact=np.zeros(len(constraints))
    for ic in range(len(constraints)):
        c_exact=np.dot(np.conj(psi),constraints[ic]['opsparse'].dot(psi)).real-constraints[ic]['cval']
        print("constraints exact=",c_exact,"qc=",c_qc[ic])


tstddev=config.getboolean('QC','tstddev')
if tstddev:
    #compute stddev
    stddev_count=int(config['QC']['stddev_count'])
    res=[]
    for i in range(stddev_count):
        v=measure_all_programs(nq,ansatz,mqcs,initial_point,backend,shots,seed+i,tsim,optimization_level)
        W_qc=measurements_to_interact(mqcs,v,pauli_interact)
        c_qc=measurements_to_constraints(mqcs,v,constraints)
        res.append([W_qc,c_qc])
    W_qc_sum=0
    c_qc_sum=np.zeros(len(constraints))
    W2_qc_sum=0
    c2_qc_sum=np.zeros(len(constraints))
    for i in range(stddev_count):
        W_qc_sum+=res[i][0]
        c_qc_sum[:]+=res[i][1][:]
        W2_qc_sum+=res[i][0]**2
        c2_qc_sum[:]+=res[i][1][:]**2
    
    W_qc_sum=W_qc_sum/stddev_count
    c_qc_sum[:]=c_qc_sum[:]/stddev_count
    W2_qc_sum=W2_qc_sum/stddev_count
    c2_qc_sum[:]=c2_qc_sum[:]/stddev_count

    print("stddev W=",math.sqrt(W2_qc_sum-W_qc_sum**2))
    for i in range(len(constraints)):
        print("stddev c[",i,"]=",math.sqrt(c2_qc_sum[i]-c_qc_sum[i]**2))


#exit()
rdmf_obj_eval=0
rdmf_cons_eval=0
tprintevals=True
printevalsevery=1

qiskit.utils.algorithm_globals.random_seed=seed
maxiter=int(config['QC']['maxiter'])    

x0=initial_point
tprintevals=False
penalty_exponent=2

if tsim and tsimcons:
    texact_expect=True
    penalty=0
    lagrange=np.zeros(len(constraints))
    method="trust-constr"
#    method="SLSQP"
    print("minimizing over parametrized qc-programs with constrained minimization with "+method+" (no derivatives, no sampling, no noise, i.e. exact expectation values)")
    eq_cons={'type': 'eq','fun' : lambda x: rdmf_cons(x)}
    if method=="trust-constr":
        res=minimize(rdmf_obj, x0, method=method, constraints=[eq_cons],tol=1e-4,options={'maxiter':10000,'verbose': 3,'disp': True,'initial_tr_radius':4*pi,'gtol':1e-4,'xtol':1e-4})
    else:
        res=minimize(rdmf_obj, x0, method=method, constraints=[eq_cons],tol=1e-4,options={'maxiter':10000,'verbose': 3,'iprint':2,'disp': True})
    print(res)
    
print("Augmented Lagrangian")
penalty=5
print("initial penalty=",penalty)

lagrange=np.zeros(len(constraints))
tprintevals=True

#augmented Lagrangian
for oiter in range(100):
    value=0
    point=[]
    nfev=0
    #unconstrainted steps
    print("Augmented Lagrangian: unconstrained subproblem")
    if tsim and not tsampling and not tnoise:
        texact_expect=True
        print("minimizing over parametrized qc-programs with augmented Lagrangian and LBFGS_B (numerical derivatives, no sampling, no noise, i.e., exact expectation values)")
        method="BFGS"
        res=minimize(rdmf_obj, x0, method=method,tol=1e-2,options={'maxiter':10000,'verbose': 2,'disp': True})
        point=res.x
        value=rdmf_obj(point)
        nfev=res.nfev
    else:
        texact_expect=False
        algo=config['QC']['algo']
        if algo=="COBYLA":
            print("minimizing over parametrized qc-programs with augmented Lagrangian and "+algo)
            optimizer = COBYLA(maxiter=maxiter,disp=True,tol=1e-2,callback=opt_callback)
            [point, value, nfev]=optimizer.optimize(num_vars=ansatz.num_parameters,objective_function=rdmf_obj,initial_point=x0)
        elif algo=="SPSA":
            print("minimizing over parametrized qc-programs with augmented Lagrangian and "+algo)
            optimizer = SPSA(maxiter=maxiter,callback=opt_callback)
            [point, value, nfev]=optimizer.optimize(num_vars=ansatz.num_parameters,objective_function=rdmf_obj,initial_point=x0)
        elif algo=="cal-SPSA":
            print("minimizing over parametrized qc-programs with augmented Lagrangian and "+algo)
            optimizer = SPSA(maxiter=maxiter,second_order=False,callback=opt_callback)
            [learning_rate,perturbation]=optimizer.calibrate(rdmf_obj,initial_point=x0,stability_constant=0, target_magnitude=None, alpha=0.602, gamma=0.101, modelspace=False)
            optimizer = SPSA(maxiter=maxiter,second_order=False,perturbation=perturbation,learning_rate=learning_rate)
            [point, value, nfev]=optimizer.optimize(num_vars=ansatz.num_parameters,objective_function=rdmf_obj,initial_point=x0)
        elif algo=="QNSPSA":
            print("minimizing over parametrized qc-programs with augmented Lagrangian and "+algo)
            fidelity = QNSPSA.get_fidelity(ansatz)
            qnspsa = QNSPSA(fidelity, maxiter=3000,callback=opt_callback)
            [point,value,nfev] = qnspsa.optimize(ansatz.num_parameters, rdmf_obj, initial_point=x0)

    print("point=",point)
    print("value=",value)
    print("nfev=",nfev)
    v=measure_all_programs(nq,ansatz,mqcs,point,backend,shots,seed,tsim,optimization_level)
    W_qc=measurements_to_interact(mqcs,v,pauli_interact)
    c_qc=measurements_to_constraints(mqcs,v,constraints)

    print("constraint violation sum(c^2)=",np.sum(c_qc**2))

    print("Augmented Lagrangian: penalty and multiplier update")
    x0=point
    #multiplier and penalty update
    for i in range(len(constraints)):
        lagrange[i]=lagrange[i]+penalty*c_qc[i]
    penalty=penalty*2
    print("new lagrange multiplier=",lagrange)
    print("new penalty=",penalty)




