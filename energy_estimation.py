from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
import Ternary_Tree as tt
import numpy as np
from qiskit import QuantumCircuit
from qiskit.algorithms.minimum_eigensolvers import NumPyMinimumEigensolver
from qiskit_nature.second_q.circuit.library import HartreeFock, UCC
from qiskit.algorithms.optimizers import L_BFGS_B, SLSQP
from qiskit.primitives import Estimator
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.second_q.algorithms import VQEUCCFactory
from qiskit_nature.second_q.circuit.library import UCC,UCCSD
from qiskit_nature.second_q.mappers import JordanWignerMapper, QubitConverter, BravyiKitaevMapper
from . import TernaryTreeMapper

from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit.quantum_info import Statevector
from qiskit.primitives.utils import (
    _circuit_key,
    _observable_key,
    bound_circuit_to_instruction,
    init_circuit,
    init_observable,
)
from qiskit.quantum_info.operators import Pauli
from qiskit.opflow import PauliTrotterEvolution, PauliOp
from .initial_state_tt import TT_initial_state

numpy_solver = NumPyMinimumEigensolver()

vqe_solver = VQEUCCFactory(Estimator(), UCC(excitations = 'sd'),  SLSQP())

def jwenergy(bond_length, at_name, nmodes =None, basis = "STO-3G"):
    driver = PySCFDriver(
        atom="H 0 0 0;" + at_name + " 0 0 " + str(bond_length),
        basis=basis,
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM,
    )
    es_problem = driver.run()
    if nmodes == None:
        nmodes = es_problem.num_spin_orbitals
        
    transformer = ActiveSpaceTransformer(num_electrons = es_problem.num_particles, 
                                         num_spatial_orbitals = nmodes//2)
    es_problem = transformer.transform(es_problem)
    
    mapper = JordanWignerMapper()
    converter = QubitConverter(mapper)
    main_operator = converter.convert(
            es_problem.second_q_ops()[0],
            num_particles=es_problem.num_particles,
            sector_locator=es_problem.symmetry_sector_locator,
        )
    qc = HartreeFock(es_problem.num_spatial_orbitals, es_problem.num_particles, converter)
    final_state = Statevector(qc)
    vqe_solver = VQEUCCFactory(Estimator(), UCC(excitations = 'sd'),  SLSQP())

    calc = GroundStateEigensolver(converter, vqe_solver)
    inenergy = final_state.expectation_value(main_operator)
    return inenergy + es_problem.nuclear_repulsion_energy

def ttenergy(bond_length, at_name, nmodes = None, basis = "STO-3G",only_hartree = False, ins = [],enum_list =None):
    
    driver = PySCFDriver(
        atom="H 0 0 0;" + at_name + " 0 0 " + str(bond_length),
#         basis="sto3g",
        basis=basis,
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM,
    )
    es_problem = driver.run()
    if nmodes == None:
        nmodes = es_problem.num_spin_orbitals
    if at_name == "Li":
        transformer = ActiveSpaceTransformer(num_electrons = 2, 
                                         num_spatial_orbitals = nmodes//2)
    else:
        transformer = ActiveSpaceTransformer(num_electrons = es_problem.num_particles, 
                                         num_spatial_orbitals = nmodes//2)
    
    es_problem = transformer.transform(es_problem)

    mapper = tt.TernaryTreeMapper(es_problem)
    converter = QubitConverter(mapper)
    main_operator = converter.convert(
            es_problem.second_q_ops()[0],
            num_particles=es_problem.num_particles,
            sector_locator=es_problem.symmetry_sector_locator,
        )
    qc = TT_initial_state(es_problem.num_spin_orbitals,es_problem.num_particles,enum_list = enum_list, ins = ins) 
    final_state = Statevector(qc)
    
    vqe_solver = VQEUCCFactory(Estimator(), UCC(excitations = 'sd'),  SLSQP(), initial_state = qc)
    calc = GroundStateEigensolver(converter, vqe_solver)
    if not only_hartree:
        res = calc.solve(es_problem)
        inenergy = final_state.expectation_value(main_operator)
    
        return res.total_energies, inenergy + res.nuclear_repulsion_energy 
    else:
        inenergy = final_state.expectation_value(main_operator)

        return inenergy + es_problem.nuclear_repulsion_energy,inenergy + es_problem.nuclear_repulsion_energy

import pyscf
def energy_classic(bond_length, at_name, nmodes = None, basis = "STO-3G"):

    mol = pyscf.M(
        atom = "H 0 0 0;" + at_name + " 0 0 " + str(bond_length),
        basis = basis)
    mf = mol.HF().run()
    mycc = mf.CISD().run()

    return mycc.e_tot, mf.e_tot

map_dict = {"JW": JordanWignerMapper, "BK": BravyiKitaevMapper, "TT": TernaryTreeMapper}

def get_depth(problem, encoding = "TT"):
    nmodes = problem.num_spin_orbitals
    if encoding == "TT":
        mapper = map_dict[encoding](problem)
        qc = TT_initial_state(problem.num_spin_orbitals, problem.num_particles) 
        converter = QubitConverter(mapper)
    else:
        mapper = map_dict[encoding]()
        converter = QubitConverter(mapper)
        qc = HartreeFock(problem.num_spatial_orbitals, problem.num_particles, converter)
    
    u = UCC(qubit_converter = converter,  num_spatial_orbitals = problem.num_spin_orbitals//2, num_particles = problem.num_particles, excitations = 'sd')

    qc = qc.compose(u.decompose().decompose().decompose())
    return qc.depth()
    
    