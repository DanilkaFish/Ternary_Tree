from __future__ import annotations

import numpy as np
from qiskit import QuantumRegister
from qiskit.circuit.library import BlueprintCircuit
from qiskit.opflow import PauliSumOp
from qiskit.utils.validation import validate_min

from qiskit_nature.second_q.mappers import BravyiKitaevSuperFastMapper, QubitConverter
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers.qubit_mapper import QubitMapper
from .pauli_tables import _pauli_table_TT, _pauli_table_JW
from qiskit.quantum_info.operators import Pauli
from qiskit.opflow import PauliTrotterEvolution, PauliOp
from numpy import pi
import copy

from .BaseTree import  st_enumeration, EnumInfo, NodeInfo
from .TernaryTree import TernaryTree
import re
dict_prod = {"II" : "I", "XX" : "I", "YY" : "I","ZZ" : "I", 
             "XY" : 'Z',"YX" : 'Z',"XZ" : 'Y',"ZX" : 'Y',"YZ" : 'X',"ZY" : 'X',
            "IX" : "X","XI" : "X","YI" : "Y","IY" : "Y","IZ" : "Z","ZI" : "Z"}

def count_non_I(syms):
    d = {'I': 0,'X': 1,'Y': 1,'Z': 1}
    n = 0
    for sym in syms:
        n += d[sym]
    return n
dict_prod_coef = {"II" : ["I",1] , "XX" :  ["I",1] , "YY" :  ["I",1] ,"ZZ" :  ["I",1] , 
             "XY" :  ["Z", 1j] ,"YX" : ["Z", -1j],"XZ" : ["Y",-1j], "ZX" : ["Y",1j],"YZ" : ["X",1j],"ZY" : ["X",-1j],
            "IX" : ["X",1], "XI" : ["X",1], "YI" : ["Y",1],"IY" : ["Y",1],"IZ" : ["Z",1],"ZI" : ["Z",1]}


def prod_exp(pt, pauli,signphi = False):
    """
    transformation pauli_table under unitary operation \exp^{i\pi/2 pauli}.
    return new puali_table with its coeff.
    """
    prodpt = []
    coef = []
    coefpt = pt[1]
    pt = pt[0]
    
    def prod(paulstr1,paulstr2,coef1 = 1):
        newpaulstr = ""
        if not signphi:
            coefpr1 = coef1*1j
        else:
            
            coefpr1 = -coef1*1j
        coefpr2 = coefpr1
        
        for i in range(len(paulstr1)):
            newpaulstr += dict_prod_coef[paulstr1[i] + paulstr2[i]][0] 
            coefpr1 = coefpr1*(dict_prod_coef[paulstr1[i] + paulstr2[i]][1])
            coefpr2 =  coefpr2*(dict_prod_coef[paulstr2[i] + paulstr1[i]][1])
        if abs((coefpr2 - coefpr1).real) < 0.0001 and abs((coefpr2 - coefpr1).imag) < 0.0001:
            return paulstr1,coef1
        else:
            return newpaulstr, coefpr2
    for i in range(len(pt)):
        newpaulstr1, coefpr1 = prod(pt[i][0],pauli,coefpt[i][0])
        newpaulstr2, coefpr2 = prod(pt[i][1],pauli,coefpt[i][1])
        prodpt.append([newpaulstr1,newpaulstr2])
        coef.append([coefpr1,coefpr2])
    return [prodpt,coef]


class InitstateTTInfo:
    def __init__(self,tt = None, nmodes = 4,num_list = None):
        """
        Class for signs check after unitary operation \exp^{i\pi/2 pauli} 
        """
        if isinstance(tt, TernaryTree):
            self.tt = tt
            

        else:
            self.tt = TernaryTree(nmodes)
        self.pttt = _pauli_table_TT(tt = tt)
        self.prod = self.prod_table(self.pttt[0])
        
    def prod_table(self, pauli_table_str):
        def prod_pauli(g1,g2):
            g1 = re.findall("[XYZI]",g1)
            g2 = re.findall("[XYZI]",g2)
            prod = [None]*len(g1)
            for i in range(len(g1)):
                prod[i] = dict_prod[g1[i] + g2[i]]
            return prod
        prod = [None]*len(pauli_table_str)
        for i, a in enumerate(pauli_table_str):
            prod[i] = ''.join(prod_pauli(pauli_table_str[i][0],pauli_table_str[i][1]) ) 

        for i, a in enumerate(prod):
            prod[i] = prod[i]
        return prod

    def check_signs(self):
        pt = self.pttt
        help_dict = {'XY': -1,'YX': 1}
        def check_order(s1,s2):
            sign = 1
            for index in range(len(s1)):

                if (s1[index] + s2[index]) in help_dict:
                    _index = index
                    sign = sign*help_dict[s1[index] + s2[index]]

                if (s1[index] + s2[index]) in ['ZI',"IZ"]:
                    if index in signs:
                        sign = -sign*signs[index]
                    else:
                        return False
            signs[_index] = int(sign*(pt[1][i][0]*pt[1][i][1]*1j).real)
            return True

        signs = {}
        i = 0
        l = len(pt[0])
        while len(signs) < l:
    #         if i not in signs:
            check_order(pt[0][i][0],pt[0][i][1])
            if i < l - 1:
                i +=1
            else:
                i = 0
        return signs
    
    def _possible_pairs(self, prod, qubit_list):
        possible_pair = []
        L = len(prod)
        for i, n in enumerate(qubit_list):
            for k in qubit_list[i + 1:]:            
                for j in range(0,L):
                    l = count_non_I(prod[j][n] + prod[j][k])
                    if l % 2 != 0:
                        break
                        
                if l % 2 == 0:
                    possible_pair.append([n,k])
        return possible_pair
   
    
    
class TT_initial_state(BlueprintCircuit):
    
    def __init__(
        self,
        nmodes = 4,
        num_particles = None,
        tt = None,
        **kwargs
    ) -> None:
        """
        Args:
            nmodes: The number of spatial orbitals.
            num_particles: The number of particles as a tuple storing the number of alpha and
                beta-spin electrons in the first and second number, respectively.
            qubit_mapper: a :class:`~qiskit_nature.second_q.mappers.QubitConverter` instance.

        Raises:
            NotImplementedError: If ``qubit_mapper`` contains
                :class:`~qiskit_nature.second_q.mappers.BravyiKitaevSuperFastMapper`. See
                https://github.com/Qiskit/qiskit-nature/issues/537 for more information.
        """

        super().__init__()

        if isinstance(tt,TernaryTree):
            self.tt = tt
        else:
            self.tt = TernaryTree(4)
        self.n_qubits = self.tt.n_qubits
        if num_particles is None:
            self.num_particles = (self.tt.nmodes//2,self.tt.nmodes//2)
        else:
            self.num_particles = num_particles
        self.qregs = [QuantumRegister(self.n_qubits, name="q")]
        self._build()
        
    def _check_configuration(self, raise_on_failure: bool = True) -> bool:
        """Check if the configuration of the HartreeFock class is valid.
        Args:
            raise_on_failure: Whether to raise on failure.
        Returns:
            True, if the configuration is valid and the circuit can be constructed. Otherwise
            returns False. Errors are only raised when raise_on_failure is set to True.

        Raises:
            ValueError: If the number of spatial orbitals is not specified or less than one.
            ValueError: If the number of particles is not specified or less than zero.
            ValueError: If the number of particles of any kind is less than zero.
            ValueError: If the number of spatial orbitals is smaller than the number of particles
                of any kind.
            ValueError: If the qubit converter is not specified.
            NotImplementedError: If the specified qubit converter is a
                :class:`~qiskit_nature.second_q.mappers.BravyiKitaevSuperFastMapper` instance.
        """
        return True
    
    def _build(self) -> None:
        """
        Construct the Hartree-Fock initial state given its parameters for ternarytree via algorithm implemented in
        TernaryTree.to0vac
        Returns:
            QuantumCircuit: a quantum circuit preparing the Hartree-Fock
            initial state given a number of spatial orbitals, number of particles and
            a qubit converter.
        """
        if self._is_built:
            return

        super()._build()

        itt = InitstateTTInfo(self.tt)
        tt = copy.deepcopy(self.tt)
        s = tt.to0vac()
        print(s)
        # Majorana fermion operators signs check after transformations
        for i in range(self.num_particles[0]):
            itt.pttt[1][i][1] = -itt.pttt[1][i][1]
        for i in range(self.num_particles[1]):
            itt.pttt[1][i + self.tt.nmodes//2][1] *= -1
        for sym in s:
            print(itt.pttt)
            itt.pttt = prod_exp(itt.pttt, sym, signphi = False)
        print(itt.pttt)
        print(tt)
        signs = itt.check_signs()
        print(signs)
#        Hartee Fock state preparation
        for k in signs:
            if signs[k] == 1:
                self.x(self.n_qubits - k - 1)
        for sym in reversed(s):
            self.compose(PauliTrotterEvolution().evolution_for_pauli(PauliOp(Pauli(sym), pi/4)).to_circuit(), inplace = True)
            
    