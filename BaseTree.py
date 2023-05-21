from __future__ import annotations
from numpy import trunc, log
import copy
import itertools

gate_name = {0: 'X', 1: 'Y', 2: 'Z'}
dict_prod = {"II" : "I", "XX" : "I", "YY" : "I","ZZ" : "I", 
             "XY" : 'Z',"YX" : 'Z',"XZ" : 'Y',"ZX" : 'Y',"YZ" : 'X',"ZY" : 'X',
            "IX" : "X","XI" : "X","YI" : "Y","IY" : "Y","IZ" : "Z","ZI" : "Z"}

class EnumInfo:
    def __init__(self,num):
        self.num = num
    def __str__(self):
        return str(self.num)
    def __repr__(self):
        return "'" + str(self.num) + "'"
    def __bool__(self):
        return False
    
    
  
class NodeInfo:
    """
    Type for nodes representation in parent_child dictionary 
    """
    def __init__(self,parent,childs = None):
        if not childs:
            self.parent = parent
            self.childs = [None]*3
        else:
            self.parent = parent
            self.childs = childs
    def __setitem__(self, key,val):
        self.childs[key] = val 
        return val
    def __getitem__(self, key):
        return self.childs[key]
    def __str__(self):
        return str([self.parent] + self.childs)
    def __repr__(self):
        return str(self.parent) + "," + str(self.childs)
    
    
def st_enumeration(nmodes):
    """
        standart numeration for alpha_beta_tree which is effecitve on my opinion 
    """
    num_list = []
    for i in range(nmodes//2):
        num_list.append([i,nmodes - i - 1])
    for i in range(nmodes + nmodes - 1, nmodes + nmodes//2 - 1, -1):
        num_list.append([i,- i + nmodes*3 - 1])
#     for i in range(nmodes):
#         num_list.append([i,2*nmodes - i - 1])
#     for i in range(nmodes + nmodes - 1, nmodes + nmodes//2 - 1, -1):
#         num_list.append([i,-i + nmodes*3 - 1])
#     num_list = [2*i for i in range(nmodes//2)]
#     num_list.append([[2*i + 1 for i in range(nmodes//2)]])
    return num_list

class BaseTree:
    """ This is class for representation ternary tree and its manipualtions"""
    def __init__(
        self,
        n_qubits: int = 0,
        parent_child: dict[int] = None,
        num_particles: int = None,
        enum_list: list[int] = None,
        nmodes: int = None,
        **kwargs
        ):
        """
        self.parent_child = {parent_number : NodeInfo} -- tree's structure
        self.enum_list = [int] -- branches numeration
        """
        
        self.parent_child = parent_child
        if parent_child:
            self.check_data()
            self.check_height()
            
        else:
            self.n_qubits = n_qubits
            if nmodes:
                self.nmodes = nmodes
            else:
                self.nmodes = n_qubits        
            self._height = 0
        
        if enum_list == None:
            self._enum_list = st_enumeration(self.nmodes)
        elif len(enum_list) == self.nmodes:
            self._enum_list = enum_list[:]
        else:
            raise ValueError
            
        l = []
        for pair in self.enum_list:
            l += pair
        self._enum_list = l

        self.max_tree = self.build_max_tree()
        
        self.max_tree = self.build_max_tree()
        if not parent_child:
            self.parent_child = copy.deepcopy(self.max_tree)
            self.build_alpha_beta_tree()
        else:
            self.set_num_qubit_transform()
            self.num_branches()
        
        
    def check_data(self):
        self.n_qubits = len(self.parent_child)
        self.nmodes = 0
        for nodes in self.parent_child:
            for child in self.parent_child[nodes].childs:
                if isinstance(child,EnumInfo): 
                    self.nmodes += 1 
        self.nmodes = self.nmodes // 2
    @property
    def enum_list(self):
        return self._enum_list
    @enum_list.setter
    def enum_list(self,num_list):
        if num_list:
            self._enum_list = num_list
        else:
            self._enum_list = st_enumeration(self.nmodes) 
        if len(self._enum_list[0]) == 2:
            l = []
            for pair in self.enum_list:
                l += pair
            self._enum_list = l
    @property
    def min_height(self):
        """
        Return minimal possible tree's height for given number of fermionic modes 
        """
        return int(trunc(log(2*self.n_qubits + 1)/log(3) - 0.001)) + 1 
    
    @property
    def height(self):
        """
        Return tree's height
        """
        self.check_height()
        return self._height
    
    def check_height(self):
        """
        This method should be wherever the height could be changed 
        """
        h = 0
        h_max = 0
        def down(parent,h):
            h+= 1
            nonlocal h_max
            for index, child in enumerate(self.parent_child[parent].childs):
                if child:
                    down(child, h)
                elif h > h_max:
                    h_max = h
        down(0,h)
        self._height = h_max
    
    @property
    def num_nodes(self):
        """
        return number of tree's nodes (number of qubits)
        """
        return len(self.parent_child)
    
    def set_num_qubit_transform(self):
        """
        function for qubit consistent numeration due to its possible deleting or inserting
        """
        self.num_qubit_transform = {}
        for i, parent in enumerate(self.parent_child):
            self.num_qubit_transform[parent] = i
            
            
    def build_alpha_beta_tree(self):
        r"""
            Supports only even number electrons. Divide tree on 2 equals parts for alpha and beta electrons.
        """
        self.parent_child = self.build_max_tree()
        num_nodes = self.num_nodes - self.n_qubits
        L = self.min_height
        l = L
        center = (3**(L-1) - 1)//2 + 3**(L-1)//2 
        parent = 0
        while (3**l - 1)//2 > num_nodes:
            parent = self.parent_child[parent][1]
            l -= 1
        left = center - 3**(l-1)//2 - 1
        right = center + 3**(l-1)//2 + 1
        
        if ((3**l - 1)//2 - num_nodes) %2 == 0:
            self.delete_node(parent)
            num_nodes = num_nodes - (3**l - 1)//2
        else:
            for child in self.parent_child[parent]:
                self.delete_node(child)
            num_nodes = num_nodes - (3**l - 1)//2 + 1
            
        i = 0
        while self.num_nodes > self.n_qubits:
            self.delete_node(right + i)
            self.delete_node(left - i)
            i += 1
        parent = 0 
        while self.parent_child[parent][1]:
            parent = self.parent_child[parent][1]
        self.parent_child[parent][1] = False
        self.set_num_qubit_transform()
        self._height = self.min_height
        self.renum_nodes()
        self.num_branches()
        
    def build_max_tree(self):
        """
        Build tree with all possible nodes and childs with height = self.min_height. Сonvenient for obtaining the necessary structures through nodes removal.
        """
        parent_child = {} 
        L = self.min_height
        full_nodes = (3**L - 1) // 2
        n = 0
        first_parent = 0
        first_child = 1
        for l in range(L - 1):
            for parent in range(first_parent, first_parent + 3**l):
                parent_child[parent] = NodeInfo(first_parent + (parent - first_parent ) // 3 - 3**(l-1))
                for child in range(3):
                    parent_child[parent][child] = first_child + child
                first_child += 3    
            first_parent = first_parent + 3**l
        for parent in range(first_parent, first_parent + 3**(L-1)):
            parent_child[parent] = NodeInfo(first_parent + (parent - first_parent ) // 3 - 3**(L-2))
        self._height = self.min_height
        return parent_child
    
    
    def renum_nodes(self):
        """
        Renumerate nodes in parent_childs after after changing the tree
        """
        self.set_num_qubit_transform()
        for parent in list(self.parent_child):
            for index, child in enumerate(self.parent_child[parent]):
                if child:
                    self.parent_child[parent][index] = self.num_qubit_transform[child]
            if self.parent_child[parent].parent >= 0:
                self.parent_child[parent].parent = self.num_qubit_transform[self.parent_child[parent].parent]
            self.parent_child[self.num_qubit_transform[parent]] = self.parent_child.pop(parent)
        self.set_num_qubit_transform()
        
        
    def num_branches(self,enum_list = None):
        """
        Numerate branches of the tree
        """
        if enum_list:
            self.enum_list = enum_list
            
        k = []
        s = []
        i = 0
        def down(parent,k):
            nonlocal s, i
            for index, child in enumerate(self.parent_child[parent].childs):
                if child:
                    down(child, k +  [[self.num_qubit_transform[parent], gate_name[index]]])
                elif child != False:
                    self.parent_child[parent][index] = EnumInfo(self.enum_list.index(i) + 1)
                    i += 1
        down(0,k)
        return self

    def delete_node(self,node):
        """
        Remove nodes and its childs from parent child
        """
        def erase(parent):
            if parent:
                for child in self.parent_child[parent]:
                    erase(child) 
                self.parent_child.pop(parent)
            
        def erase_from_parent(child):
            parent = self.parent_child[child].parent
            if parent in self.parent_child:
                i = self.parent_child[parent].childs.index(child)
                self.parent_child[parent].childs[i] = None
        if node in self.parent_child:
            erase_from_parent(node)
            erase(node)
            
    
    
    
    def branches(self, get_num = False):
        """
        convert parent_child to list of branches
        """
        k = []
        s = []
        num_list = []
        def down(parent,k):
            nonlocal s
            for index, child in enumerate(self.parent_child[parent].childs):
                if child:
                    down(child, k +  [[self.num_qubit_transform[parent], gate_name[index]]])
                elif child != False:
                    s.append(k + [[self.num_qubit_transform[parent], gate_name[index]]])
                    num_list.append(child.num)
        down(0,k)
        if get_num:
            return s, num_list
        else:
            return s
    
    
    def branch_transposition(self,first_node,first_edge, second_node, second_edge):
        """
        Transpose subtree or branches  
        """
        s = ["I"]*len(self.parent_child)
        closest_node = self.closest_parent(first_node,second_node)
        node1 = first_node
        edge1 = gate_name[first_edge]
        while closest_node != node1:
            s[node1] = edge1 
            for index,child in enumerate(self.parent_child[self.parent_child[node1].parent]):
                if child == node1:
                    edge1 = gate_name[index]
            node1 = self.parent_child[node1].parent
        node2 = second_node
        edge2 = gate_name[second_edge]
        while closest_node != node2:
            s[node2] = edge2 
            for index,child in enumerate(self.parent_child[self.parent_child[node2].parent]):
                if child == node2:
                    edge2 = gate_name[index]
            node2 = self.parent_child[node2].parent
        s[closest_node] = dict_prod[edge1 + edge2]
        
        # меняю информацию у first_node: child := parent_shild[second_node][second_edge] и [second_node][second_edge] если нужно
        info = self.parent_child[first_node][first_edge]
        self.parent_child[first_node][first_edge] = self.parent_child[second_node][second_edge] 
        if self.parent_child[second_node][second_edge]: # здесь еще нужно поменять инфо у [second_node][second_edge]
            old_child = self.parent_child[second_node][second_edge]
            self.parent_child[old_child].parent = first_node

        # меняю информацию у second_node: child := parent_shild[first_node][first_edge] и [first_node][first_edge] если нужно
        self.parent_child[second_node][second_edge] = info
        if info:
            old_child = info
            self.parent_child[old_child].parent = second_node        

        self.check_height()
        
        return ''.join(s)
    
    
    def closest_parent(self,first_node,second_node):
        def check_node(init_node,search_node):
            flag = True
            def down(node,search_node):
                nonlocal flag
                for index, child in enumerate(self.parent_child[node].childs):
                    if child == search_node:
                        flag = False
                    if child:
                        down(child, search_node)
            if init_node ==second_node:
                return False
            down(init_node,search_node)           
            return flag
        while check_node(first_node,second_node):
            first_node = self.parent_child[first_node].parent
        return first_node
    
    def __str__(self):
        k = []
        s = []
        L = self.height
        enum_list = []
        def down(parent,k):
            nonlocal s, enum_list
            for index, child in enumerate(self.parent_child[parent].childs):
                if child:
                    down(child, k +  [str(parent) + gate_name[index]])
                elif child != False:
                    s.append( k + [str(parent) +  gate_name[index]] )
                    enum_list.append(self.parent_child[parent][index])
        
        down(0,k)
        pr = ''
        num_space = 0
        k = 4
        for l in range(0,L):
            for branch in s:
                if len(branch) > l:
                    pr += branch[l]+ ' '*(k-len(branch[l]))
                else:
                    pr += " "*k
            pr = pr + '\n'
        for num in enum_list:
            pr += str(num) + ' '*(k-len(str(num)))
                                  
        return pr + '\n ---------------------------------------------------------------'
    
    
#     old too difficult output

#     def __str__(self):
#         k = []
#         s = []
#         def down(parent,k):
#             nonlocal s
#             for index, child in enumerate(self.max_tree[parent].childs):
#                 if child:
#                     down(child, k +  [str(parent) + gate_name[index]])
#                 else:
#                     s.append( k + [str(parent) +  gate_name[index]] )
#         down(0,k)
#         L = self.height
#         pr = ''
#         num_space = 0

#         for l in range(L-1,-1, -1):
#             num_space = 0
#             space = '_'*num_space
#             h = ''
#             sym = '_'
#             for index, branch in enumerate(s):
#                 i = index # индекс для вывода номера ветви
#                 h += space
#                 if (index )% 3 == 0:
#                     h += " "
#                 if len(branch) > l and index%3**(L - l - 1) == 3**(L-l - 1)//2 and int(branch[l][:-1]) in self.parent_child:
#                     h += branch[l][-1]
#                 else:
                    
#                     h += ' '
#                 num_space = num_space * 3 + 2
#             pr = "\n" + h + pr
#         return pr
    