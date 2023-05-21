# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""A product formula base for decomposing non-commuting operator exponentials."""

from typing import Callable, Optional, Union, Any, Dict
from functools import partial
import numpy as np
from qiskit.circuit.parameterexpression import ParameterExpression
from qiskit.circuit.quantumcircuit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Pauli

from copy import copy

from qiskit_nature.second_q.circuit.library import HartreeFock, UCC
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

def optimize_ucc(
    u: QuantumCircuit
):
    n_qubits = u.num_qubits
    ins = [inst for inst in u]
    pass_manager = generate_preset_pass_manager(optimization_level = 3)
#     print(u)
#     print(u.decompose(reps = 2))
    pauli_list = [inst.operation.name[7:-1] for inst in ins]
#     print(pauli_list)
    p,targ, r = pauli_composing(pauli_list)
    u0 = QuantumCircuit(n_qubits)
#     print(pauli_list)
    for i,pauli in enumerate(p):
        param = ins[pauli_list.index(pauli)].operation.params[0]
        qc = evolve_pauli(
            pauli = pauli,
            time = param,
            cx_structure = "fountain",
            label = None,
            target = targ[i]
        )
        u0 = u0.compose(qc)
#     print(u0.num_parameters)
#     print(u0)
    u0  = u0.decompose(reps = 4)
    u0 = pass_manager.run(u0)   
    u0  = u0.decompose(reps = 4)
#     print(u0)
    return u0,r 

def compare(pauli1,pauli2,pos_tar):
    count = 0
    flag = False
    possible_target = []
    
    target = - 1
    if pos_tar > 0:
        if pauli1[pos_tar] == "I" or pauli2[pos_tar] == "I":
            return 0 , -1
        else: 
            for i,gate in enumerate(pauli1):
                if (gate == pauli2[i]) and not(gate == "I") and pos_tar != i:
                    count += 2
            return count, pos_tar
    else:       
        for i,gate in enumerate(pauli1):
            if (gate == pauli2[i]) and not(gate == "I"):
                count += 2
                if target < 0:
                    target = i
                    count -= 2
            elif gate + pauli2[i] in ["XY","YX"]:
                if target > 0:
                    count += 2
                target = i


#     if target == -1:
#         print("-1")
    return count, target

def pauli_composing(pl):
#     создаю копию pl
    pauli_list = [pauli for pauli in pl]
    pauli_reorder = []
    targ = [-1]
    pauli_reorder.append(pauli_list[0])
    pauli_list.pop(0)
    r = 0
    while len(pauli_list) > 0:
        index = 0
        count = 0
        target = -1
        for i, pauli in enumerate(pauli_list):
            count0, target0 = compare(pauli_reorder[-1],pauli,targ[-1])
            if  count0 > count:
                count = count0
                index = i
                target = target0

        r+=count
        pauli_reorder.append(pauli_list[index])
        i = 1
        targ.append(target)
        if target >= 0: 
            while i < len(targ) and targ[-1 - i] < 0:
                targ[-1 - i] = target
                i += 1
        pauli_list.pop(index)

    return [pauli_reorder, targ, r]

def count_length(pauli_list):
    l = 0 
    for pauli in pauli_list:
        l += len(pauli)
        for gate in pauli:
            if gate == "I":
                l-=1
        l -= 1
    return l*2


r"""Pauli exponentioation redefining"""
def evolve_pauli(
    pauli: str,
    time: Union[float, ParameterExpression] = 1.0,
    cx_structure: str = "chain",
    label: Optional[str] = None,
    target = None
) -> QuantumCircuit:
    r"""Construct a circuit implementing the time evolution of a single Pauli string.

    For a Pauli string :math:`P = \{I, X, Y, Z\}^{\otimes n}` on :math:`n` qubits and an
    evolution time :math:`t`, the returned circuit implements the unitary operation

    .. math::

        U(t) = e^{-itP}.

    Since only a single Pauli string is evolved the circuit decomposition is exact.

    Args:
        pauli: The Pauli to evolve.
        time: The evolution time.
        cx_structure: Determine the structure of CX gates, can be either "chain" for
            next-neighbor connections or "fountain" to connect directly to the top qubit.
        label: A label for the gate.

    Returns:
        A quantum circuit implementing the time evolution of the Pauli.
    """
    if target == -1:
        target = None
    pauli = Pauli(pauli)
    num_non_identity = len([label for label in pauli.to_label() if label != "I"])
    
    # first check, if the Pauli is only the identity, in which case the evolution only
    # adds a global phase
    if num_non_identity == 0:
        definition = QuantumCircuit(pauli.num_qubits, global_phase=-time)
    # if we evolve on a single qubit, if yes use the corresponding qubit rotation
    elif num_non_identity == 1:
        definition = _single_qubit_evolution(pauli, time)
    # same for two qubits, use Qiskit's native rotations
    elif num_non_identity == 2:
        definition = _two_qubit_evolution(pauli, time, cx_structure, target)
    # otherwise do basis transformation and CX chains
    else:
        definition = _multi_qubit_evolution(pauli, time, cx_structure, target)

    definition.name = f"exp(it {pauli.to_label()})"

    return definition


def _single_qubit_evolution(pauli, time):
    definition = QuantumCircuit(pauli.num_qubits)
    # Note that all phases are removed from the pauli label and are only in the coefficients.
    # That's because the operators we evolved have all been translated to a SparsePauliOp.
    for i, pauli_i in enumerate(reversed(pauli.to_label())):
        if pauli_i == "X":
            definition.rx(2 * time, i)
        elif pauli_i == "Y":
            definition.ry(2 * time, i)
        elif pauli_i == "Z":
            definition.rz(2 * time, i)

    return definition


def _two_qubit_evolution(pauli, time, cx_structure, target = None):
    # Get the Paulis and the qubits they act on.
    # Note that all phases are removed from the pauli label and are only in the coefficients.
    # That's because the operators we evolved have all been translated to a SparsePauliOp.
    labels_as_array = np.array(list(reversed(pauli.to_label())))
    qubits = np.where(labels_as_array != "I")[0]
    if qubits[0] != target:
        qubits[0],qubits[1] = qubits[1],qubits[0]
    labels = np.array([labels_as_array[idx] for idx in qubits])

    definition = QuantumCircuit(pauli.num_qubits)

    # go through all cases we have implemented in Qiskit
#     if all(labels == "X"):  # RXX
#         definition.rxx(2 * time, qubits[0], qubits[1])
#     elif all(labels == "Y"):  # RYY
#         definition.ryy(2 * time, qubits[0], qubits[1])
#     elif all(labels == "Z"):  # RZZ
#         definition.rzz(2 * time, qubits[0], qubits[1])
#     elif labels[0] == "Z" and labels[1] == "X":  # RZX
#         definition.rzx(2 * time, qubits[0], qubits[1])
#     elif labels[0] == "X" and labels[1] == "Z":  # RXZ
#         definition.rzx(2 * time, qubits[1], qubits[0])
#     else:  # all the others are not native in Qiskit, so use default the decomposition
    definition = _multi_qubit_evolution(pauli, time, cx_structure, target)
    return definition


def _multi_qubit_evolution(pauli, time, cx_structure, target = None):
    # get diagonalizing clifford
    cliff = diagonalizing_clifford(pauli)

    # get CX chain to reduce the evolution to the top qubit
    if cx_structure == "chain":
        chain = cnot_chain(pauli)
    else:
        chain = cnot_fountain(pauli,target)
    
    # determine qubit to do the rotation on
#     target = None
    # Note that all phases are removed from the pauli label and are only in the coefficients.
    # That's because the operators we evolved have all been translated to a SparsePauliOp.
    if target == None:
        for i, pauli_i in enumerate(reversed(pauli.to_label())):
            if pauli_i != "I":
                target = i
                break

    # build the evolution as: diagonalization, reduction, 1q evolution, followed by inverses
    definition = QuantumCircuit(pauli.num_qubits)
    definition.compose(cliff, inplace=True)
    definition.compose(chain, inplace=True)
    definition.rz(2 * time, target)
    definition.compose(chain.inverse(), inplace=True)
    definition.compose(cliff.inverse(), inplace=True)

    return definition


def diagonalizing_clifford(pauli: Pauli) -> QuantumCircuit:
    """Get the clifford circuit to diagonalize the Pauli operator.

    Args:
        pauli: The Pauli to diagonalize.

    Returns:
        A circuit to diagonalize.
    """
    cliff = QuantumCircuit(pauli.num_qubits)
    for i, pauli_i in enumerate(reversed(pauli.to_label())):
        if pauli_i == "Y":
            cliff.sdg(i)
        if pauli_i in ["X", "Y"]:
            cliff.h(i)

    return cliff


def cnot_chain(pauli: Pauli) -> QuantumCircuit:
    """CX chain.

    For example, for the Pauli with the label 'XYZIX'.

    .. parsed-literal::

                       ┌───┐
        q_0: ──────────┤ X ├
                       └─┬─┘
        q_1: ────────────┼──
                  ┌───┐  │
        q_2: ─────┤ X ├──■──
             ┌───┐└─┬─┘
        q_3: ┤ X ├──■───────
             └─┬─┘
        q_4: ──■────────────

    Args:
        pauli: The Pauli for which to construct the CX chain.

    Returns:
        A circuit implementing the CX chain.
    """

    chain = QuantumCircuit(pauli.num_qubits)
    control, target = None, None

    # iterate over the Pauli's and add CNOTs
    for i, pauli_i in enumerate(pauli.to_label()):
        i = pauli.num_qubits - i - 1
        if pauli_i != "I":
            if control is None:
                control = i
            else:
                target = i

        if control is not None and target is not None:
            chain.cx(control, target)
            control = i
            target = None

    return chain


def cnot_fountain(pauli: Pauli, target = None) -> QuantumCircuit:
    """CX chain in the fountain shape.

    For example, for the Pauli with the label 'XYZIX'.

    .. parsed-literal::

             ┌───┐┌───┐┌───┐
        q_0: ┤ X ├┤ X ├┤ X ├
             └─┬─┘└─┬─┘└─┬─┘
        q_1: ──┼────┼────┼──
               │    │    │
        q_2: ──■────┼────┼──
                    │    │
        q_3: ───────■────┼──
                         │
        q_4: ────────────■──

    Args:
        pauli: The Pauli for which to construct the CX chain.

    Returns:
        A circuit implementing the CX chain.
    """
    if target == None:
        chain = QuantumCircuit(pauli.num_qubits)
        control, target = None, None
        for i, pauli_i in enumerate(reversed(pauli.to_label())):
            if pauli_i != "I":
                if target is None:
                    target = i
                else:
                    control = i

            if control is not None and target is not None:
                chain.cx(control, target)
                control = None
    else:
        
        chain = QuantumCircuit(pauli.num_qubits)
        control = None
        for i, pauli_i in enumerate(reversed(pauli.to_label())):
            if pauli_i != "I":
                if i != target:
                    control = i

            if control is not None and target is not None:
                chain.cx(control, target)
                control = None

    return chain


def _default_atomic_evolution(operator, time, cx_structure):
    if isinstance(operator, Pauli):
        # single Pauli operator: just exponentiate it
        evolution_circuit = evolve_pauli(operator, time, cx_structure)
    else:
        # sum of Pauli operators: exponentiate each term (this assumes they commute)
        pauli_list = [(Pauli(op), np.real(coeff)) for op, coeff in operator.to_list()]
        name = f"exp(it {[pauli.to_label() for pauli, _ in pauli_list]})"
        evolution_circuit = QuantumCircuit(operator.num_qubits, name=name)
        for pauli, coeff in pauli_list:
            evolution_circuit.compose(evolve_pauli(pauli, coeff * time, cx_structure), inplace=True)

    return evolution_circuit
