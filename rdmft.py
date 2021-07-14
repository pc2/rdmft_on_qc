import os
try:
    import numpy as np
except ImportError:
    print("installing numpy...")
    os.system("pip3 install numpy")
    print("installed numpy. Please restart.")
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
    from qiskit.algorithms import VQE
    from qiskit.algorithms.optimizers import L_BFGS_B,SPSA
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
    from qiskit import Aer
    from qiskit_nature.algorithms import GroundStateEigensolver
    from qiskit_nature.converters.second_quantization import QubitConverter
    from qiskit_nature.drivers import PySCFDriver, UnitsType, Molecule
    from qiskit_nature.mappers.second_quantization import ParityMapper
    from qiskit_nature.problems.second_quantization import ElectronicStructureProblem
    from qiskit_nature.problems.second_quantization.electronic import ElectronicStructureProblem
    from qiskit_nature.converters.second_quantization import QubitConverter
    from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper    
    
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


#VQE in qiskit from https://github.com/Qiskit/qiskit-nature/
np.random.seed(999999)

# Use PySCF, a classical computational chemistry software
# package, to compute the one-body and two-body integrals in
# electronic-orbital basis, necessary to form the Fermionic operator

molecule = Molecule(geometry=[['H', [0., 0., 0.]],['H', [0., 0., 0.735]]],charge=0, multiplicity=1)

driver = PySCFDriver(molecule=molecule,unit=UnitsType.ANGSTROM,basis='sto3g')
problem = ElectronicStructureProblem(driver)

# generate the second-quantized operators
second_q_ops = problem.second_q_ops()
print("second_q_ops=",second_q_ops)

main_op = second_q_ops[0]

num_particles = (problem.molecule_data_transformed.num_alpha,problem.molecule_data_transformed.num_beta)
num_spin_orbitals = 2 * problem.molecule_data.num_molecular_orbitals

#qubit_converter = QubitConverter(mapper=ParityMapper(), two_qubit_reduction=True)
qubit_converter = QubitConverter(mapper=JordanWignerMapper())
qubit_op = qubit_converter.convert(main_op,num_particles=num_particles)
print("qubit_op=",qubit_op)

# setup the initial state for the ansatz
#entangler_map = [(0, 1), (1, 2), (2, 0)]
ansatz = TwoLocal(4,rotation_blocks = ['rx', 'ry'], entanglement_blocks = 'cz',entanglement='linear', reps=2, parameter_prefix = 'y')

print(ansatz)

ansatz.draw(output='text',filename="ansatz.txt")
ansatz.draw(output='mpl',filename="ansatz.png")


# set the backend for the quantum computation
backend = Aer.get_backend('aer_simulator_statevector')

#optimizer = L_BFGS_B()
optimizer = SPSA(maxiter=100,callback=opt_callback)
print(optimizer.print_options())
initial_point = np.random.random(ansatz.num_parameters)

# setup and run VQE
algorithm = VQE(ansatz=ansatz,optimizer=optimizer,quantum_instance=backend,initial_point=initial_point)

calc = GroundStateEigensolver(qubit_converter, algorithm)
res = calc.solve(problem)

print(res)


#build 1-RDM qubit op

#build ee-interaction qubit op

#read in density matrix

#augmented Lagrangian
