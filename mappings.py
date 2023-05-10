from .pauli_tables import pauli_table_JW, pauli_table_BK, pauli_table_TT
from .pauli_weight import pauli_single_count, pauli_double_count
import copy
import numpy as np
from qiskit.quantum_info.operators import Pauli

map_dict = {"JW": pauli_table_JW, "BK": pauli_table_BK, "TT": pauli_table_TT}
pauli_table_JW
def getinf(nmodes, num_alpha, num_beta, mapper_name = "JW", flag_transposition = False,**kwargs):
    if flag_transposition:
        list_branches = get_all_branches_from_pauli_tables(nmodes, mapper_name)
    else:
        list_branches = [puali_table_to_branches(map_dict[mapper_name](nmodes, **kwargs))]
    return supinf(list_branches, num_alpha, num_beta)
    
def get_all_branches_from_pauli_tables(nmodes, mapper_name = "JW"):
    pauli_table = map_dict[mapper_name](nmodes) 
    list_tables = enumerate_pauli_table(pauli_table)
    for i in range(len(list_tables)):
        list_tables[i] = puali_table_to_branches(list_tables[i])
    return list_tables

def supinf(all_branches, num_alpha, num_beta):
    from math import factorial

    sup1 = 0
    tree_sup1 = []
    inf1 = 1000
    tree_inf1 = []
    num_inf1 = []
    
    sup2 = 0
    tree_sup2 = []
    inf2 = 1000
    tree_inf2 = []
    num_inf2 = []
    
    sup = 0
    tree_sup = []
    inf = 1000
    tree_inf = []

    for branches in all_branches:
        s1 = pauli_single_count(branches, num_alpha, num_beta)
        s2 = pauli_double_count(branches, num_alpha, num_beta)

        s = (s1[0] + s2[0])/(s1[1] + s2[1])
        s1 = s1[0]/s1[1]
        s2 = s2[0]/s2[1]
        if sup < s:
            sup = s
            tree_sup = [branches]
        if sup == s:
            tree_sup.append(branches)

        if inf > s:
            inf = s
            tree_inf = [branches]
        if inf == s:
            tree_inf.append(branches)

        if sup2 < s2:
            sup2 = s2
            tree_sup2 = [branches]
        if sup2 == s2:
            tree_sup2.append(branches)

        if inf2 > s2:
            inf2 = s2
            tree_inf2 = [branches]
        if inf2 == s2:
            tree_inf2.append(branches)
            
        if sup1 < s1:
            sup1 = s1
            tree_sup2 = [branches]
        if sup1 == s1:
            tree_sup1.append(branches)

        if inf1 > s1:
            inf1 = s1
            tree_inf1 = [branches]
        if inf1 == s1:
            tree_inf1.append(branches)
    return [inf,inf1,inf2]

# def enumerate_pauli_table(pauli_table):
#     print(len(pauli_table))
    
#     s = []
#     k = []
#     lim = 6
#     def numerate(pauli_table1,k):
#         flag = True
#         for i, branch in enumerate(pauli_table1):
#             pauli_table0 = copy.copy(pauli_table1)
#             pauli_table0.pop(i)
#             numerate(pauli_table0, k + [pauli_table1[i]])
#             flag = False
#         if flag:
#             for i in range(lim,len(pauli_table)):
#                 k += [pauli_table[i]]
#             s.append(k)
            
        
#     numerate(pauli_table[:lim], k)
# #     print(s[10])
#     return s
    
def enumerate_pauli_table(pauli_table):
    print(len(pauli_table))
    
    s = []
    k1 = []
    k2 = []
    lim = min(len(pauli_table)//2, 6)
    pauli_table_half1 = pauli_table[:len(pauli_table)//2]     
    pauli_table_half2 = pauli_table[len(pauli_table)//2:]   
    def numerate(pauli_table1, pauli_table2, k1,k2):
        flag = True
        for i, branch in enumerate(pauli_table1):
            pauli_table10 = copy.copy(pauli_table1)
            pauli_table10.pop(i)
            pauli_table20 = copy.copy(pauli_table2)
            pauli_table20.pop(i)
            numerate(pauli_table10, pauli_table20, k1 + [pauli_table1[i]],k2 + [pauli_table2[i]])
            flag = False
            
        if flag:
            for i in range(lim,len(pauli_table)//2):
                k1 += [pauli_table_half1[i]]
                k2 += [pauli_table_half2[i]]
#             s.append(k)
            s.append(k1 + k2)
#     print(pauli_table_half1)
#     print(pauli_table_half2)
    numerate(pauli_table_half1[:lim],pauli_table_half2[:lim], k1,k2)
#     print(s[0])
    return s

def puali_table_to_branches(pauli_table):
    n_qubits = len(pauli_table)
    branches = []
    for a in pauli_table:
        g1 = []
        g2 = []
        for n in range(n_qubits):
            if  a[0][n].to_label() != 'I':
                g1.append([n, a[0][n].to_label()])
                g2.append([n, a[1][n].to_label()])
        branches.append(g1)
        branches.append(g2)
    return branches
