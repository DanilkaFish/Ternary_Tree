from __future__ import annotations
from numpy import trunc, log
import copy
import itertools

gate_name = {0: 'X', 1: 'Y', 2: 'Z'}
dict_prod = {"II" : "I", "XX" : "I", "YY" : "I","ZZ" : "I", 
             "XY" : 'Z',"YX" : 'Z',"XZ" : 'Y',"ZX" : 'Y',"YZ" : 'X',"ZY" : 'X',
            "IX" : "X","XI" : "X","YI" : "Y","IY" : "Y","IZ" : "Z","ZI" : "Z"}

class EnumInfo:
    """
        Type for numeration of the branches. This is special type for tree leaves.
    """
    def __init__(self,num):
        self.num = num
    def __str__(self):
        return str(self.num)
    def __repr__(self):
        return "'" + str(self.num) + "'"
    def __bool__(self):
        return False
    def __eq__(self,other):
        if isinstance(other,EnumInfo):
            return self.num == other.num
        else:
            return False
    
class NodeInfo:
    """
    Type for nodes representation in parent_child dictionary 
    """
    def __init__(self,parent = 0,childs = None):
        """
        self.parent = int -- number of parent node
        self.childs = list(int | EnumInfo | None) -- int - number child node, EnumInfo - number of the branch, None - if there is no child
        """
        if isinstance(parent,(int,float)):
            self.parent = parent
        else:
            self.parent = 0
        if not childs:
            self._childs = [None]*3
        else:
            self._childs = childs
            
            
    @property
    def childs(self):
        return self._childs
    
    @childs.setter
    def childs(self):
        for i, child in enumerate(self.childs):
            if not isinstance(child, (int,EnumInfo)):
                self.childs[i] = False
        while len(self.childs)<3:
            self.childs.append(False)

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
        Standart numeration for alpha_beta_tree which is effecitve on my opinion
    """
    num_list = []
    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(nmodes//2):
            num_list.append([i,nmodes - i - 1])
        for i in range(nmodes + nmodes - 1, nmodes + nmodes//2 - 1, -1):
            num_list.append([i,- i + nmodes*3 - 1])
    else:
        num_list = [[2*i,2*i + 1] for i in range(nmodes)]
#     for i in range(nmodes):
#         num_list.append([i,2*nmodes - i - 1])
#     for i in range(nmodes + nmodes - 1, nmodes + nmodes//2 - 1, -1):
#         num_list.append([i,-i + nmodes*3 - 1])
#     num_list = [2*i for i in range(nmodes//2)]
#     num_list.append([[2*i + 1 for i in range(nmodes//2)]])
    return num_list



class BaseTree:
    """ This is the class for representation ternary tree and its manipualtions"""
    def __init__(
        self,
        n_qubits: int = 0,
        parent_child: dict[int] = None,
        num_particles: int = None,
        nmodes: int = None,
        enum_list: list[int] = None,
        **kwargs
        ):
        """
        self.n_qubits = int -- number of qubits or number of tree nodes
        self.parent_child = {parent_number : NodeInfo} -- tree's structure
        self.num_paritcles = tuple(int,int) -- number alpha and beta electrons respectively (Need only for branches' numeration)
        self.nmodes = int number -- numbet fermionic modes (Need only for branches' numeration)
        self.enum_list = [int] -- branches numeration
        """
        
        self._parent_child = copy.deepcopy(parent_child)
        if self._parent_child:
            self.parent_child = self._parent_child
        else:
            self.n_qubits = n_qubits
            if nmodes:
                self.nmodes = nmodes
            else:
                self.nmodes = n_qubits        
            self._height = 0
            self.build_alpha_beta_tree()
        self.enum_list = enum_list        
    
    @property
    def parent_child(self):
        return self._parent_child
    
    @parent_child.setter
    def parent_child(self,pch):
        self._check_parent_child(pch)
        self._parent_child = pch
        self.set_data()
        self.check_height()
                
    @property
    def enum_list(self):
        """
        Branches numeration
        """
        return self._enum_list
    
    @enum_list.setter
    def enum_list(self,num_list):
        if num_list:
            self._enum_list = num_list
        else:
            self._enum_list = st_enumeration(self.nmodes) 
        if isinstance(self._enum_list[0], list) :
            l = []
            for pair in self.enum_list:
                l += pair
            self._enum_list = l
        self.num_branches()
        
    @property
    def min_height(self):
        """
        Return minimal possible tree's height for given number of qubits 
        """
        return int(trunc(log(2*self.n_qubits + 1)/log(3) - 0.001)) + 1 
    
    @property
    def height(self):
        """
        Return tree's height
        """
        self.check_height()
        return self._height
    
    def _check_parent_child(self,parent_child):
        "Small check of parent_child object"
#             parent_child = self.parent_child
        check = {}
        edges = {}
        def append(parent):
            if not isinstance(parent, int):
                raise ValueError("Parent should be int object, not " + str(type(parent)))
            if parent in check:
                raise TreeStructureError('Multiple initialization of parent node with number ' + str(parent))
            check[parent] = None 

        for parent in parent_child:
            append(parent)
            if isinstance(parent_child[parent], NodeInfo):
                pass
#                 for child in self.parent_childs[parent].childs:
#                     if isinstance(child, int):
#                         l = []
#                         edges[parent] = edges.get(parent, l).append(child)
            else:
                raise ValueError("KeyValue in self.parent_child should be NodeInfo object, not " + str(type( parent_child[parent])))
           
    def set_data(self):
        """
        Set n_qubits and nmodes from parent_child info
        """
        self.n_qubits = len(self.parent_child)
        self.nmodes = 0
        num_list = []
        renum_flag = False
        for nodes in self.parent_child:
            for child in self.parent_child[nodes].childs:
                if isinstance(child,EnumInfo): 
                    self.nmodes += 1
                    num_list.append(child.num)
                if child is None:
                    renum_flag = True
                    self.nmodes += 1
                        
        if self.nmodes % 2 == 1:
            raise TreeStructureError("parent_child should contain only even number branches")
        self.nmodes = self.nmodes // 2
        if renum_flag or any([not num in range(0,self.nmodes*2) for num in num_list]):
            self.enum_list = st_enumeration(self.nmodes)
        else:
            self.enum_list = num_list
            
    def check_height(self):
        """
        This method should be wherever the height could be changed, Set BaseTree._height
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
    
    
    def set_num_qubit_transform(self):
        """
        Function for qubit consistent numeration due to its possible deleting or inserting.
        """
        self.num_qubit_transform = {}
        for i, parent in enumerate(self.parent_child):
            self.num_qubit_transform[parent] = i
            
    def renum_nodes(self):
        """
        Renumerate nodes in parent_childs after after changing the tree
        """
        self.set_num_qubit_transform()
        for parent in list(self.parent_child):
            # Changing KeyValue
            for index, child in enumerate(self.parent_child[parent].childs):
                if child:
                    self.parent_child[parent][index] = self.num_qubit_transform[child]
            if self.parent_child[parent].parent >= 0:
                self.parent_child[parent].parent = self.num_qubit_transform[self.parent_child[parent].parent]
            # Creating new renumerated key and deleting old key
            self.parent_child[self.num_qubit_transform[parent]] = self.parent_child.pop(parent)
        self.set_num_qubit_transform()
        
    # Old method for tree inizialization       
    def build_alpha_beta_tree(self, height = 0):
        r"""
            Supports only even number electrons. Divide tree on 2 equals parts for alpha and beta electrons.
        """
        def _delete_node(self,node):
            """
            Remove nodes and its childs from parent child
            """
            def erase(parent):
                if parent:
                    for child in parent_child[parent].childs:
                        if child:
                            erase(child) 
                    parent_child.pop(parent)

            def erase_from_parent(child):
                parent = parent_child[child].parent
                if parent in parent_child:
                    i = parent_child[parent].childs.index(child)
                    parent_child[parent].childs[i] = None
            if node in parent_child:
                erase_from_parent(node)
                erase(node)
        
        parent_child = self.build_max_tree(height)
        num_nodes = len(parent_child) - self.n_qubits
        L = self.min_height
        l = L
        center = (3**(L-1) - 1)//2 + 3**(L-1)//2 
        parent = 0
        while (3**l - 1)//2 > num_nodes:
            parent = parent_child[parent][1]
            l -= 1
        left = center - 3**(l-1)//2 - 1
        right = center + 3**(l-1)//2 + 1
        
        if ((3**l - 1)//2 - num_nodes) %2 == 0:
            _delete_node(self,parent)
            num_nodes = num_nodes - (3**l - 1)//2
        else:
            for child in parent_child[parent]:
                _delete_node(self,child)
            num_nodes = num_nodes - (3**l - 1)//2 + 1
            
        i = 0
        while len(parent_child) > self.n_qubits:
            _delete_node(self,right + i)
            _delete_node(self,left - i)
            i += 1
        parent = 0 
        while parent_child[parent][1]:
            parent = parent_child[parent][1]
        parent_child[parent][1] = False
        self.parent_child = parent_child
        self.set_num_qubit_transform()
        self.renum_nodes()
        
    def build_max_tree(self, height = 0):
        """
        Build tree with all possible nodes and childs with height = self.min_height | height. Сonvenient for obtaining the necessary structures through nodes removal.
        """
        parent_child = {} 
        L = max(height,self.min_height)
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
        return parent_child
    
        
    def num_branches(self,enum_list = None):
        """
        Consistently numerate branches of the tree
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
                    down(child, k +  [[parent, gate_name[index]]])
                elif child != False:
                    self.parent_child[parent][index] = EnumInfo(self.enum_list.index(i) + 1)
                    i += 1
        down(0,k)
        return self

    def delete_node(self,node):
        """
        Remove nodes and its childs from parent child
        """
        flag = False
        def erase(parent):
            nonlocal flag
            if parent:
                for child in self.parent_child[parent].childs:
                    if child:
                        erase(child)
                    elif child == False:
                        flag = True
                self.parent_child.pop(parent)
            
        def parent_ch_num(child):
            parent = self.parent_child[child].parent
            if parent in self.parent_child:
                i = self.parent_child[parent].childs.index(child)
                self.parent_child[parent].childs[i] = None
                return parent, i
        if node in self.parent_child:
            parent, i = parent_ch_num(node)
            erase(node)
            if flag:
                self.parent_child[parent].childs[i] = False
                
        self.renum_nodes()
        self.set_data()
        self.check_height()
        self.num_branches()
        print(self)
    
    
    def branches(self, get_num = False):
        """
        Convert parent_child to list of branches. Used by pauli tables
        """
        k = []
        s = []
        num_list = []
        def down(parent,k):
            nonlocal s
            for index, child in enumerate(self.parent_child[parent].childs):
                if child:
                    down(child, k +  [[parent, gate_name[index]]])
                elif child != False:
                    s.append(k + [[parent, gate_name[index]]])
                    num_list.append(child.num)
        down(0,k)
        if get_num:
            return s, num_list
        else:
            return s    
    
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
    
   
    