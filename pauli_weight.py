from __future__ import annotations
from numpy import trunc, log
import copy
import itertools

def pauli_single_count(tree_branches, num_alpha, num_beta):
    ngamma = len(tree_branches)
    a_alpha_oc = [tuple([tree_branches[i], tree_branches[i + 1]]) for i in range(num_alpha)]
    a_beta_oc = [tuple([tree_branches[i + ngamma//2], tree_branches[i + 1 + ngamma//2]]) for i in range(num_beta)]
    a_alpha_un = [tuple([tree_branches[i], tree_branches[i + 1]]) for i in range(num_alpha,ngamma//4)]
    a_beta_un = [tuple([tree_branches[i + ngamma//2], tree_branches[i + 1 + ngamma//2]]) for i in range(num_beta,ngamma//4)]
    def num_in_prod(branch1,branch2):
        num = len(branch1) + len(branch2)
        for el1 in branch1:
            for el2 in branch2:
                if el1[0] == el2[0]:
                    if el1[1] == el2[1]:
                        num -= 1
                    num -= 1
        return num
    sum_weight = 0
    n = 0
    for a_alpha_oc_el in a_alpha_oc:
        for a_alpha_un_el in a_alpha_un:
            n += 1
            sum_weight += num_in_prod(a_alpha_oc_el[0],a_alpha_un_el[0])
            n += 1
            sum_weight += num_in_prod(a_alpha_oc_el[1],a_alpha_un_el[1])
    for a_beta_oc_el in a_beta_oc:
        for a_beta_un_el in a_beta_un:
            n += 1
            sum_weight += num_in_prod(a_beta_oc_el[0],a_beta_un_el[0])
            n += 1
            sum_weight += num_in_prod(a_beta_oc_el[1],a_beta_un_el[1])
    return sum_weight, n


def pauli_double_count(tree_branches, num_alpha, num_beta):
    ngamma = len(tree_branches)
    a_alpha_oc = [tuple([tree_branches[i], tree_branches[i + 1]]) for i in range(num_alpha)]
    a_beta_oc = [tuple([tree_branches[i + ngamma//2], tree_branches[i + 1 + ngamma//2]]) for i in range(num_beta)]
    a_alpha_un = [tuple([tree_branches[i], tree_branches[i + 1]]) for i in range(num_alpha,ngamma//4)]
    a_beta_un = [tuple([tree_branches[i + ngamma//2], tree_branches[i + 1 + ngamma//2]]) for i in range(num_beta,ngamma//4)]
    
    def num_in_prod(branches):
        num = sum([len(branch) for branch in branches])
        for n in range(ngamma//2):
            mas = []
            for branch in branches:
                i = 0
                l = len(branch)
                while (i < l) and (branch[i][0]<= n):
                    if branch[i][0] == n:
                        mas.append(branch[i])
                    i += 1
            if len(mas) == 2:
                if mas[0][1] == mas[1][1]:
                    num -= 2
                else:
                    num -= 1
            if len(mas) == 3:
                num -=3
                if (mas[0][1] == mas[1][1]) or (mas[2][1] == mas[1][1]) or (mas[0][1] == mas[2][1]):
                    num += 1
            
            if len(mas) == 4:
                num -=3
                if (mas[0][1] == mas[1][1] and mas[2][1] == mas[3][1]) or (mas[0][1] == mas[2][1] 
                    and mas[1][1] == mas[3][1]) or (mas[0][1] == mas[3][1] and mas[1][1] == mas[2][1]):
                    num -= 1
                
        return num
    sum_weight = 0
    n = 0
    for a_alpha_oc_el1, a_alpha_oc_el2 in itertools.product(a_alpha_oc,a_alpha_oc):
        for a_alpha_un_el1,a_alpha_un_el2 in itertools.product(a_alpha_un,a_alpha_un):
            mas = [a_alpha_oc_el1, a_alpha_oc_el2, a_alpha_un_el1, a_alpha_un_el2]
            if not((a_alpha_oc_el1 == a_alpha_oc_el2) or (a_alpha_un_el1 == a_alpha_un_el2)):
                mas0 = [i[0] for i in mas]
                for i in range(4):
                    mas0[i] = mas[i][1]
                    sum_weight += num_in_prod(mas0)
                n += 4
                mas0 = [i[1] for i in mas]
                for i in range(4):
                    mas0[i] = mas[i][0]
                    sum_weight += num_in_prod(mas0)
                n += 4
                
    for a_alpha_oc_el1, a_alpha_oc_el2 in itertools.product(a_alpha_oc,a_beta_oc):
        for a_alpha_un_el1,a_alpha_un_el2 in itertools.product(a_alpha_un,a_beta_un):
            mas = [a_alpha_oc_el1, a_alpha_oc_el2, a_alpha_un_el1, a_alpha_un_el2]
            if not((a_alpha_oc_el1 == a_alpha_oc_el2) or (a_alpha_un_el1 == a_alpha_un_el2)):
                mas0 = [i[0] for i in mas]
                for i in range(4):
                    mas0[i] = mas[i][1]
                    sum_weight += num_in_prod(mas0)
                n += 4
                mas0 = [i[1] for i in mas]
                for i in range(4):
                    mas0[i] = mas[i][0]
                    sum_weight += num_in_prod(mas0)
                n += 4
                
    for a_alpha_oc_el1, a_alpha_oc_el2 in itertools.product(a_beta_oc,a_beta_oc):
        for a_alpha_un_el1,a_alpha_un_el2 in itertools.product(a_beta_un,a_beta_un):
            mas = [a_alpha_oc_el1, a_alpha_oc_el2, a_alpha_un_el1, a_alpha_un_el2]
            if not((a_alpha_oc_el1 == a_alpha_oc_el2) or (a_alpha_un_el1 == a_alpha_un_el2)):
                mas0 = [i[0] for i in mas]
                for i in range(4):
                    mas0[i] = mas[i][1]
                    sum_weight += num_in_prod(mas0)
                n += 4
                mas0 = [i[1] for i in mas]
                for i in range(4):
                    mas0[i] = mas[i][0]
                    sum_weight += num_in_prod(mas0)
                n += 4
                
    return sum_weight, n 
