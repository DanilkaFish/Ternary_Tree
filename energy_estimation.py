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
from . import TernaryTreeMapper, EnumInfo, NodeInfo, TernaryTree
from .my_evolution import optimize_ucc
from qiskit.algorithms.minimum_eigensolvers import VQE
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
from .parent_childs import get_bin_parent_child
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
def ttenergy_opt4( bond_length = "1",at_name = "H",nmodes = 4, basis="STO-3G", active_orbitals = None, only_hartree = True):
    if not active_orbitals:
        active_orbitals = [i for i in range(nmodes//2)]
    driver = PySCFDriver(
        atom="H 0 0 0;" + at_name + " 0 0 " + str(bond_length),
    #         basis="sto3g",
        basis=basis,
        charge=0,
        spin=0,
        unit=DistanceUnit.ANGSTROM,
    )
    es_problem = driver.run()
    transformer = ActiveSpaceTransformer(num_electrons = es_problem.num_alpha + es_problem.num_beta , 
                                         num_spatial_orbitals = nmodes//2  , active_orbitals = active_orbitals)
    es_problem = transformer.transform(es_problem)

    parent_child = get_bin_parent_child(nmodes)
    tt = TernaryTree(parent_child = parent_child)
#     print(tt)
#     print(tt.parent_child)
    mapper = TernaryTreeMapper(es_problem, tt = tt)
    converter = QubitConverter(mapper)
    # print(converter)
    qc = TT_initial_state(es_problem.num_spin_orbitals, es_problem.num_particles, tt = tt)
    u = UCC(qubit_converter = converter,  num_spatial_orbitals = es_problem.num_spin_orbitals//2, 
        num_particles = (es_problem.num_alpha, es_problem.num_beta), excitations = 'sd')
    u = u.decompose(reps = 2)
    # print(u)

    ucc, _ =optimize_ucc(u) 
    
    converter = QubitConverter(mapper)
    main_operator = converter.convert(
            es_problem.second_q_ops()[0],
            num_particles=es_problem.num_particles,
            sector_locator=es_problem.symmetry_sector_locator,
        )
    final_state = Statevector(qc)
    
    
    
    if not only_hartree:
        ucc = qc.compose(ucc)
        vqe_solver = VQE(Estimator(), ucc,  SLSQP())
        calc = GroundStateEigensolver(converter, vqe_solver)
        res = calc.solve(es_problem)
        inenergy = final_state.expectation_value(main_operator)
    
        return res.total_energies, inenergy + res.nuclear_repulsion_energy 
    else:
        inenergy = final_state.expectation_value(main_operator)

        return inenergy + es_problem.nuclear_repulsion_energy,inenergy + es_problem.nuclear_repulsion_energy
    
    
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

    mapper = tt.TernaryTreeMapper(es_problem, enum_list = enum_list)
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
    
    