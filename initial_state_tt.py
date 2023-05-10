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
             "XY" :  ["Z",1j] ,"YX" : ["Z", - 1j],"XZ" : ["Y",-1j], "ZX" : ["Y",1j],"YZ" : ["X",1j],"ZY" : ["X",-1j],
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
    def __init__(self,nmodes,num_list = None):
        """
        Class for signs check after unitary operation \exp^{i\pi/2 pauli} 
        """
        self.nmodes = nmodes 
        self.num_list = num_list
        self.oqe, self.tqe = self.get_entenglement() #[tuple(1,"X")] , [tuple([n,m], ["IX","IY"...])]
        self.pttt = _pauli_table_TT(self.nmodes, self.num_list)
        self.ptjw = _pauli_table_JW(self.nmodes, self.num_list)
        self.prod = self.prod_table(_pauli_table_TT(nmodes,self.num_list)[0])
        
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

#     def check_signs(self):
#         pt = self.pttt
#         def check_order(s):
#             for index,sym in enumerate(s):
#                 if sym == "X":
#                     return [1,index]
#                 if sym == "Y":
#                     return [-1,index]
#         signs = []
#         for i in range(len(pt[0])):
#             h_sign, index = check_order(pt[0][i][0])
#             sign = int((pt[1][i][0]*1j*pt[1][i][1]).real) * h_sign
#             signs.append([sign,self.nmodes - index - 1])
#         return signs
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

    def get_entenglement(self):
        prod = self.prod_table(_pauli_table_TT(self.nmodes, self.num_list)[0])
        qubit_list = [i for i in range(self.nmodes)]
        oqe = self._one_qubit_entenglement(prod, qubit_list)
        tqe = self._two_qubits_entanglement(qubit_list, prod)
        return oqe, tqe
        
    def _one_qubit_entenglement(self, prod, qubit_list):
        n_qubits = len(qubit_list)
        L = len(prod)
        sym_qub = []
        for n in qubit_list:
            sym = ''
            for j in range(L):
                if prod[j][n] not in sym:
                    sym += prod[j][n]
            if len(sym) < 3:
                sym_qub.append(tuple([n ,re.findall("[^I]",sym)]))
        for n in sym_qub:
            qubit_list.remove(n[0])
        return  sym_qub
    
    def _two_qubits_entanglement(self, qubit_list, prod):
        possible_pairs = copy.deepcopy( self._possible_pairs(prod,qubit_list) )
        tqe = []
        L = 0
        def _get_dict():
            diff_number_pair = {}
            nonlocal L 
            L = len(possible_pairs)
            
            for j in range(L):
                pair = possible_pairs[j]
                diff_number_pair[pair[0]] = diff_number_pair.get(pair[0],[]) + [pair]
                diff_number_pair[pair[1]] = diff_number_pair.get(pair[1],[]) + [pair]
            return diff_number_pair
        
        diff_number_pair = _get_dict()
        while len(diff_number_pair):
            min1, key1, min2, key2 = [self.nmodes]*4

            for key in diff_number_pair:
                l = len(diff_number_pair[key])
                if l < min1:
                    min1 = l
                    key1 = key
            for _pair in diff_number_pair[key1]:
                i = _pair.index(key1)
                l = len(diff_number_pair[_pair[i - 1]])
                if l <min2:
                    min2 = l
                    key2 = _pair[i - 1]
                    pair = _pair
            tqe.append(tuple([pair, self.get_str_pairs(prod,pair)]))
            for i in range(L-1,-1,-1):
                if pair[0] in possible_pairs[i] or pair[1] in possible_pairs[i]:
                    possible_pairs.pop(i)
            diff_number_pair = _get_dict()
        return tqe
    
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
    
    def get_str_pairs(self,prod,pair):
        str_pair = []
        L = len(prod)
        sym_qubs = []
        for j in range(L):
            sym_qubs.append(prod[j][pair[0]] + prod[j][pair[1]])
        return sym_qubs
    
    


    
class TT_initial_state(BlueprintCircuit):
    
    def __init__(
        self,
        nmodes,
        num_particles,
        enum_list = None,
        ins = None,
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
        self.enum_list = enum_list
        self.ins = ins
        self.nmodes = nmodes
        self.num_particles = num_particles
        self.qregs = [QuantumRegister(nmodes, name="q")]
        print(self.enum_list)
        self.tt = kwargs.get("tt", TernaryTree(self.nmodes, enum_list = self.enum_list))
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
        Construct the Hartree-Fock initial state given its parameters.
        Returns:
            QuantumCircuit: a quantum circuit preparing the Hartree-Fock
            initial state given a number of spatial orbitals, number of particles and
            a qubit converter.
        """
        if self._is_built:
            return

        super()._build()
#         tt = TernaryTree(self.nmodes)
#         s = tt.tojw()
#         itt = InitstateTTInfo(self.nmodes)
#         for sym in s:
#             itt.pttt = prod_exp(itt.pttt, sym, signphi = False)
#         signs = itt.check_signs()
#         for i in range(self.num_particles[0]):
#             signs[i][0] *= - 1
#         for i in range(self.num_particles[1]):
#             signs[i + self.nmodes//2][0] *= - 1
#         for k in signs:
#             if k[0] == -1:
#                 self.x(k[1])
#         for sym in reversed(s):
#             self.compose(PauliTrotterEvolution().evolution_for_pauli(PauliOp(Pauli(sym), pi/4)).to_circuit(), inplace = True)
        
        s = self.tt.to0vac()
        itt = InitstateTTInfo(self.nmodes)
#       
        for i in range(self.num_particles[0]):
            itt.pttt[1][i][1] = -itt.pttt[1][i][1]
        for i in range(self.num_particles[1]):
            itt.pttt[1][i + self.nmodes//2][1] *= -1
        for sym in s:
            itt.pttt = prod_exp(itt.pttt, sym, signphi = False)
        
        
        signs = itt.check_signs()
#         Подготовка состояния Хартри-Фока
        for k in signs:
            if signs[k] == 1:
                self.x(self.nmodes - k - 1)
        
        for sym in reversed(s):
#             sym1 = ''
#             for r in reversed(sym):
#                 sym1 += r
            self.compose(PauliTrotterEvolution().evolution_for_pauli(PauliOp(Pauli(sym), pi/4)).to_circuit(), inplace = True)
            
    