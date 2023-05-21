from __future__ import annotations

from functools import lru_cache


from numpy import trunc, log
from qiskit.opflow import PauliSumOp
from qiskit.quantum_info.operators import Pauli, SparsePauliOp

from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers import FermionicMapper
import copy
import numpy as np
from .TernaryTree import TernaryTree
from qiskit_nature import QiskitNatureError
from qiskit_nature.second_q.operators import SparseLabelOp
from .pauli_tables import pauli_table_TT

class TernaryTreeMapper(FermionicMapper):  # pylint: disable=missing-class-docstring
    
    def __init__(self, es_problem,trans = None,enum_list = None, tt = 0):
        """The Teranry Tree fermion-to-qubit mapping. Standart mapping is alpha_beta_tree"""
#         print(type(tt),type(TernaryTree))
        if isinstance(tt,TernaryTree):
        
            self.tt = tt
            self.tt.num_branches()
            self.nmodes = tt.nmodes
            if trans:
                self.trans = trans
            else:
                self.trans = [i for i in range(self.nmodes)]
        self.num_alpha, self.num_beta = es_problem.num_particles
        self.n_electrons = self.num_alpha + self.num_beta
        
        super().__init__(allows_two_qubit_reduction=False)

    @property
    def tree_height(self):
        return int(trunc(log(2*self.nmodes + 1)/log(3) - 0.001)) + 1

    def pauli_table(self, nmodes):
        pt = copy.copy(pauli_table_TT(nmodes, tt = self.tt))
        pauli_table = [None]*self.nmodes 
        for i,t in enumerate(self.trans):
            pauli_table[t] = pt[i]
        
        return pauli_table
    
    @lru_cache(maxsize=32)
    def sparse_pauli_operators(cls, nmodes: int) -> tuple[list[SparsePauliOp], list[SparsePauliOp]]:
        """Generates the cached :class:`.SparsePauliOp` terms.

        This uses :meth:`.QubitMapper.pauli_table` to construct a list of operators used to
        translate the second-quantization symbols into qubit operators.

        Args:
            nmodes: the number of modes for which to generate the operators.

        Returns:
            Two lists stored in a tuple, consisting of the creation and annihilation  operators,
            applied on the individual modes.
        """
        times_creation_op = []
        times_annihilation_op = []

        for paulis in cls.pauli_table(nmodes):
            real_part = SparsePauliOp(paulis[0], coeffs=[0.5])
            imag_part = SparsePauliOp(paulis[1], coeffs=[0.5j])

            # The creation operator is given by 0.5*(X - 1j*Y)
            creation_op = real_part - imag_part
            times_creation_op.append(creation_op)

            # The annihilation operator is given by 0.5*(X + 1j*Y)
            annihilation_op = real_part + imag_part
            times_annihilation_op.append(annihilation_op)
            n_qubits = len(paulis[0])
        
        return times_creation_op, times_annihilation_op, n_qubits
    
#     @classmethod
    def mode_based_mapping(cls, second_q_op: SparseLabelOp, nmodes: int) -> PauliSumOp:
        """Utility method to map a `SparseLabelOp` to a `PauliSumOp` using a pauli table.

        Args:
            second_q_op: the `SparseLabelOp` to be mapped.
            nmodes: the number of modes for which to generate the operators.

        Returns:
            The `PauliSumOp` corresponding to the problem-Hamiltonian in the qubit space.

        Raises:
            QiskitNatureError: If number length of pauli table does not match the number
                of operator modes, or if the operator has unexpected label content
        """
#         print("nmodes = ", nmodes)
        times_creation_op, times_annihilation_op, n_qubits = cls.sparse_pauli_operators(nmodes)
#         print("times_creation_op = ", times_creation_op, "\n", times_annihilation_op)
        # make sure ret_op_list is not empty by including a zero op
        # Here are differencies
        ret_op_list = [SparsePauliOp("I" * n_qubits, coeffs=[0])]

        for terms, coeff in second_q_op.terms():
            # 1. Initialize an operator list with the identity scaled by the `coeff`
            ret_op = SparsePauliOp("I" * n_qubits, coeffs=np.array([coeff]))
#             print("terms = ",terms)
            # Go through the label and replace the fermion operators by their qubit-equivalent, then
            # save the respective Pauli string in the pauli_str list.
            for term in terms:
                
                char = term[0]
                if char == "":
                    break
                position = int(term[1])
                if char == "+":
                    ret_op = ret_op.compose(times_creation_op[position], front=True)
                    
                elif char == "-":
                    ret_op = ret_op.compose(times_annihilation_op[position], front=True)
                # catch any disallowed labels
                else:
                    raise QiskitNatureError(
                        f"FermionicOp label included '{char}'. Allowed characters: I, N, E, +, -"
                    )
#                 print("ret_op = ", ret_op)
            ret_op_list.append(ret_op)

        return PauliSumOp(SparsePauliOp.sum(ret_op_list).simplify())

    def map(self, second_q_op: FermionicOp) -> PauliSumOp:
        return self.mode_based_mapping(second_q_op, second_q_op.register_length)

    
class JordanWignerMapper(FermionicMapper):  # pylint: disable=missing-class-docstring
    def __init__(self,trans = None):
        """The Jordan-Wigner fermion-to-qubit mapping."""
        self.trans = trans
        super().__init__(allows_two_qubit_reduction=False)

#     @classmethod
#     @lru_cache(maxsize=32)
    def pauli_table(self, nmodes: int) -> list[tuple[Pauli, Pauli]]:
        pauli_table = []

        for i in range(nmodes):
            a_z = np.asarray([1] * i + [0] + [0] * (nmodes - i - 1), dtype=bool)
            a_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
            b_z = np.asarray([1] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
            b_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
            # c_z = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
            # c_x = np.asarray([0] * nmodes, dtype=bool)
            pauli_table.append((Pauli((a_z, a_x)), Pauli((b_z, b_x))))
            # TODO add Pauli 3-tuple to lookup table
        pt = copy.deepcopy(pauli_table)
        if self.trans:
            for i,t in enumerate(self.trans):
                pauli_table[t] = pt[i]
        return pauli_table
    
    def sparse_pauli_operators(cls, nmodes: int) -> tuple[list[SparsePauliOp], list[SparsePauliOp]]:
        """Generates the cached :class:`.SparsePauliOp` terms.

        This uses :meth:`.QubitMapper.pauli_table` to construct a list of operators used to
        translate the second-quantization symbols into qubit operators.

        Args:
            nmodes: the number of modes for which to generate the operators.

        Returns:
            Two lists stored in a tuple, consisting of the creation and annihilation  operators,
            applied on the individual modes.
        """
        times_creation_op = []
        times_annihilation_op = []

        for paulis in cls.pauli_table(nmodes):
            real_part = SparsePauliOp(paulis[0], coeffs=[0.5])
            imag_part = SparsePauliOp(paulis[1], coeffs=[0.5j])

            # The creation operator is given by 0.5*(X - 1j*Y)
            creation_op = real_part - imag_part
            times_creation_op.append(creation_op)

            # The annihilation operator is given by 0.5*(X + 1j*Y)
            annihilation_op = real_part + imag_part
            times_annihilation_op.append(annihilation_op)

        return (times_creation_op, times_annihilation_op)
    
    def mode_based_mapping(cls, second_q_op: SparseLabelOp, nmodes: int) -> PauliSumOp:
        """Utility method to map a `SparseLabelOp` to a `PauliSumOp` using a pauli table.

        Args:
            second_q_op: the `SparseLabelOp` to be mapped.
            nmodes: the number of modes for which to generate the operators.

        Returns:
            The `PauliSumOp` corresponding to the problem-Hamiltonian in the qubit space.

        Raises:
            QiskitNatureError: If number length of pauli table does not match the number
                of operator modes, or if the operator has unexpected label content
        """
        times_creation_op, times_annihilation_op = cls.sparse_pauli_operators(nmodes)

        # make sure ret_op_list is not empty by including a zero op
        ret_op_list = [SparsePauliOp("I" * nmodes, coeffs=[0])]

        for terms, coeff in second_q_op.terms():
            # 1. Initialize an operator list with the identity scaled by the `coeff`
            ret_op = SparsePauliOp("I" * nmodes, coeffs=np.array([coeff]))

            # Go through the label and replace the fermion operators by their qubit-equivalent, then
            # save the respective Pauli string in the pauli_str list.
            for term in terms:
                char = term[0]
                if char == "":
                    break
                position = int(term[1])
                if char == "+":
                    ret_op = ret_op.compose(times_creation_op[position], front=True)
                elif char == "-":
                    ret_op = ret_op.compose(times_annihilation_op[position], front=True)
                # catch any disallowed labels
                else:
                    raise QiskitNatureError(
                        f"FermionicOp label included '{char}'. Allowed characters: I, N, E, +, -"
                    )
            ret_op_list.append(ret_op)

        return PauliSumOp(SparsePauliOp.sum(ret_op_list).simplify())

    def map(self, second_q_op: FermionicOp) -> PauliSumOp:
        return self.mode_based_mapping(second_q_op, second_q_op.register_length)
    
    
class BravyiKitaevMapper(FermionicMapper):  # pylint: disable=missing-class-docstring
    def __init__(self, trans = None):
        """The Bravyi-Kitaev fermion-to-qubit mapping."""
        self.trans = trans
        super().__init__(allows_two_qubit_reduction=False)
    def sparse_pauli_operators(cls, nmodes: int) -> tuple[list[SparsePauliOp], list[SparsePauliOp]]:
        """Generates the cached :class:`.SparsePauliOp` terms.

        This uses :meth:`.QubitMapper.pauli_table` to construct a list of operators used to
        translate the second-quantization symbols into qubit operators.

        Args:
            nmodes: the number of modes for which to generate the operators.

        Returns:
            Two lists stored in a tuple, consisting of the creation and annihilation  operators,
            applied on the individual modes.
        """
        times_creation_op = []
        times_annihilation_op = []

        for paulis in cls.pauli_table(nmodes):
            real_part = SparsePauliOp(paulis[0], coeffs=[0.5])
            imag_part = SparsePauliOp(paulis[1], coeffs=[0.5j])

            # The creation operator is given by 0.5*(X - 1j*Y)
            creation_op = real_part - imag_part
            times_creation_op.append(creation_op)

            # The annihilation operator is given by 0.5*(X + 1j*Y)
            annihilation_op = real_part + imag_part
            times_annihilation_op.append(annihilation_op)

        return (times_creation_op, times_annihilation_op)
    
    def mode_based_mapping(cls, second_q_op: SparseLabelOp, nmodes: int) -> PauliSumOp:
        """Utility method to map a `SparseLabelOp` to a `PauliSumOp` using a pauli table.

        Args:
            second_q_op: the `SparseLabelOp` to be mapped.
            nmodes: the number of modes for which to generate the operators.

        Returns:
            The `PauliSumOp` corresponding to the problem-Hamiltonian in the qubit space.

        Raises:
            QiskitNatureError: If number length of pauli table does not match the number
                of operator modes, or if the operator has unexpected label content
        """
        times_creation_op, times_annihilation_op = cls.sparse_pauli_operators(nmodes)

        # make sure ret_op_list is not empty by including a zero op
        ret_op_list = [SparsePauliOp("I" * nmodes, coeffs=[0])]

        for terms, coeff in second_q_op.terms():
            # 1. Initialize an operator list with the identity scaled by the `coeff`
            ret_op = SparsePauliOp("I" * nmodes, coeffs=np.array([coeff]))

            # Go through the label and replace the fermion operators by their qubit-equivalent, then
            # save the respective Pauli string in the pauli_str list.
            for term in terms:
                char = term[0]
                if char == "":
                    break
                position = int(term[1])
                if char == "+":
                    ret_op = ret_op.compose(times_creation_op[position], front=True)
                elif char == "-":
                    ret_op = ret_op.compose(times_annihilation_op[position], front=True)
                # catch any disallowed labels
                else:
                    raise QiskitNatureError(
                        f"FermionicOp label included '{char}'. Allowed characters: I, N, E, +, -"
                    )
            ret_op_list.append(ret_op)

        return PauliSumOp(SparsePauliOp.sum(ret_op_list).simplify())
    
    def pauli_table(cls, nmodes: int) -> list[tuple[Pauli, Pauli]]:
        def parity_set(j, n):
            """
            Computes the parity set of the j-th orbital in n modes.

            Args:
                j (int) : the orbital index
                n (int) : the total number of modes

            Returns:
                numpy.ndarray: Array of mode indices
            """
            indices = np.array([])
            if n % 2 != 0:
                return indices

            if j < n / 2:
                indices = np.append(indices, parity_set(j, n / 2))
            else:
                indices = np.append(
                    indices, np.append(parity_set(j - n / 2, n / 2) + n / 2, n / 2 - 1)
                )
            return indices

        def update_set(j, n):
            """
            Computes the update set of the j-th orbital in n modes.

            Args:
                j (int) : the orbital index
                n (int) : the total number of modes

            Returns:
                numpy.ndarray: Array of mode indices
            """
            indices = np.array([])
            if n % 2 != 0:
                return indices
            if j < n / 2:
                indices = np.append(indices, np.append(n - 1, update_set(j, n / 2)))
            else:
                indices = np.append(indices, update_set(j - n / 2, n / 2) + n / 2)
            return indices

        def flip_set(j, n):
            """
            Computes the flip set of the j-th orbital in n modes.

            Args:
                j (int) : the orbital index
                n (int) : the total number of modes

            Returns:
                numpy.ndarray: Array of mode indices
            """
            indices = np.array([])
            if n % 2 != 0:
                return indices
            if j < n / 2:
                indices = np.append(indices, flip_set(j, n / 2))
            elif n / 2 <= j < n - 1:
                indices = np.append(indices, flip_set(j - n / 2, n / 2) + n / 2)
            else:
                indices = np.append(
                    np.append(indices, flip_set(j - n / 2, n / 2) + n / 2), n / 2 - 1
                )
            return indices

        pauli_table = []
        # FIND BINARY SUPERSET SIZE
        bin_sup = 1
        while nmodes > np.power(2, bin_sup):
            bin_sup += 1
        # DEFINE INDEX SETS FOR EVERY FERMIONIC MODE
        update_sets = []
        update_pauli = []

        parity_sets = []
        parity_pauli = []

        flip_sets = []

        remainder_sets = []
        remainder_pauli = []
        for j in range(nmodes):

            update_sets.append(update_set(j, np.power(2, bin_sup)))
            update_sets[j] = update_sets[j][update_sets[j] < nmodes]

            parity_sets.append(parity_set(j, np.power(2, bin_sup)))
            parity_sets[j] = parity_sets[j][parity_sets[j] < nmodes]

            flip_sets.append(flip_set(j, np.power(2, bin_sup)))
            flip_sets[j] = flip_sets[j][flip_sets[j] < nmodes]

            remainder_sets.append(np.setdiff1d(parity_sets[j], flip_sets[j]))

            update_pauli.append(Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool))))
            parity_pauli.append(Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool))))
            remainder_pauli.append(
                Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
            )
            for k in range(nmodes):
                if np.in1d(k, update_sets[j]):
                    update_pauli[j].x[k] = True
                if np.in1d(k, parity_sets[j]):
                    parity_pauli[j].z[k] = True
                if np.in1d(k, remainder_sets[j]):
                    remainder_pauli[j].z[k] = True

            x_j = Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
            x_j.x[j] = True
            y_j = Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
            y_j.z[j] = True
            y_j.x[j] = True
            pauli_table.append(
                (
                    parity_pauli[j] & x_j & update_pauli[j],
                    remainder_pauli[j] & y_j & update_pauli[j],
                )
            )

        # PauliList has the phase information.
        # Here, phase is unnecessary, so the following removes phase.
        for pauli1, pauli2 in pauli_table:
            pauli1.phase = 0
            pauli2.phase = 0
        pt = copy.deepcopy(pauli_table)
        if cls.trans:
            for i,t in enumerate(cls.trans):
                pauli_table[t] = pt[i]
        return pauli_table

    def map(self, second_q_op: FermionicOp) -> PauliSumOp:
        return self.mode_based_mapping(second_q_op, second_q_op.register_length)