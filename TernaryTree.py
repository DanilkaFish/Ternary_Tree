from .BaseTree import BaseTree, EnumInfo


gate_name = {0: 'X', 1: 'Y', 2: 'Z'}
gate_index = {'X': 0,'Y': 1, 'Z': 2}
dict_prod = {"II" : "I", "XX" : "I", "YY" : "I","ZZ" : "I", 
             "XY" : 'Z',"YX" : 'Z',"XZ" : 'Y',"ZX" : 'Y',"YZ" : 'X',"ZY" : 'X',
            "IX" : "X","XI" : "X","YI" : "Y","IY" : "Y","IZ" : "Z","ZI" : "Z"}

class TernaryTree(BaseTree):
    """
    Ternary tree object for initial state preparation
    """
        
    def branch_transposition(self,first_node,first_edge, second_node, second_edge):
        """
        Transpose subtrees or branches  wtih preserving of nodes numeration. Return pauli operator implementing this transformation.
        first_node = int -- node number in tree
        first_edge = 0|1|2 -- edge which connected first_node with its parent 
        second_node = int -- node number in tree
        second_edge = 0|1|2 -- edge which connected second_node with its parent
        """
        s = ["I"]*len(self.parent_child)
        closest_node, _ = self._closest_parent(first_node,second_node)
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
        self.set_data()
        return ''.join(s)
    
    
    def _closest_parent(self,first_node,second_node):
        """
        Return closest parent to the first and second node
        """
        def check_node(init_node,search_node):
            """
            Check whether init_node is parent of search_node or not
            """
            flag = True
            dist = 0
            def down(node,search_node,_dist):
                nonlocal flag, dist
                for index, child in enumerate(self.parent_child[node].childs):
                    if child == search_node:
                        dist = _dist + 1
                        flag = False
                    if child:
                        down(child, search_node,_dist + 1)
            if init_node == second_node:
                return False, 0
            down(init_node,search_node,0)           
            return flag, dist
        flag, dist0 = check_node(first_node,second_node) 
        i = 0
        while flag:
            first_node = self.parent_child[first_node].parent
            flag, dist0 = check_node(first_node,second_node) 
            i += 1
        return first_node, i + dist0
    
    
    def to0vac(self):
        
        """
        Use algorithm described in the graduate work
        """
        help_dict = {0:1, 1:0}
        def _find_second_branch(node,edge):
            """
            This function search branch to tranpose with given 
            """
            parent2,flag = -1, False
            def down(parent):
                nonlocal parent2, flag
                for index, child in enumerate(self.parent_child[parent].childs):
                    if flag:
                        break
                    if parent == node and index == edge:
                        i=0
                    elif isinstance(child,EnumInfo) and index == 2:                        
                        parent2 = parent
                        flag = True
                    elif child:
                        down(child)
            down(0)
            return parent2

        def _check_xy():
            """
            find unoccupied branch for transposition with Z...Z branch  
            """
            parentf = 0 
            indexf = 0
            def down(parent):
                nonlocal parentf,indexf 
                for index, child in enumerate(self.parent_child[parent].childs):
                    if isinstance(child, bool):
                        parentf, indexf = parent, index 
                    if child:
                        down(child)
            down(0)
            return [parentf,indexf]
        
        def check_pair(pair):
            if min(pair[0],pair[1]) % 2 == 1 and abs(pair[0] - pair[1]) == 1:
                return True
            else:
                return False
        def find_trans(pair_list,index):
            pair = pair_list[index]
            indexes = []
            if check_pair(pair):
                return False
            for i, ipair in enumerate(pair_list[index + 1:]):
                if check_pair([ipair[0],pair[0]]):
                    return [2*index + 1, 2*(i + index + 1) ]
                if check_pair([ipair[0],pair[1]]):
                    return [2*index, 2*(i + index + 1) ]
                if check_pair([ipair[1],pair[0]]):
                    return [2*index + 1, 2*(i + index + 1) + 1]
                if check_pair([ipair[1],pair[1]]):
                    return [2*index, 2*(i + index + 1) + 1]
       
                
        def get_pair(self,node1,edge1, node2, edge2, branches):
            l = [branch[-1] for branch in branches]
            index = self.enum_list[l.index([node1,gate_name[edge1]])]
            if index % 2 ==0:
                _index = index + 1
            else:
                _index = index - 1
            _node, _edge = branches[self.enum_list.index(_index)][-1]
            _, dist = self._closest_parent(node2,_node) 
            return _node, gate_index[_edge], dist
        
        
        
        s = []
        parent = 0 
        
        # find Z...Z branch to eliminate it
        child = self.parent_child[parent][2]
        while child:
            parent = child
            child = self.parent_child[parent][2]
        if not isinstance(child,bool):
            first, index = _check_xy()
            s.append(self.branch_transposition(first, index, parent,2))
        
#         Лучше через parent_child сделать
#         first leaf
        i = 0
        for i in range(self.nmodes):
            branches = self.branches()
            node, edge = branches[self.enum_list.index(2*i)][-1]
            edge = gate_index[edge]
            if edge == 0:
                if isinstance(self.parent_child[node][1],bool):
                    s.append(self.branch_transposition(node,1,node,2))
                if isinstance(self.parent_child[node][1],EnumInfo):
                    _node = node
                    _edge = 1
                
                else:
                    _node = self.parent_child[node][1]
                    _edge = 2
                    while not isinstance(self.parent_child[_node][2],EnumInfo):
                        if isinstance(self.parent_child[_node][2],bool):
                            s.append(self.branch_transposition(_node,2,_node,1))
                        else:
                            _node = self.parent_child[_node][2]
                branches = self.branches()
                
                node1,edge1,dist1 = get_pair(self,node, edge, _node, _edge, branches)
#                 node2,edge2,dist2 = get_pair(self,_node, _edge, node, edge, branches) 
                if node1 != _node or edge1 != _edge:
                    
#                     if dist1 < dist2:
                    s.append(self.branch_transposition(_node,_edge,node1,edge1))
#                     else:
#                         s.append(self.branch_transposition(node,edge,node2,edge2))
            if edge == 1:
                if isinstance(self.parent_child[node][0],bool):
                    s.append(self.branch_transposition(_node,0,node,2))
                if isinstance(self.parent_child[node][0],EnumInfo):
                    _node = node
                    _edge = 0
                
                else:
                    _node = self.parent_child[node][0]
                    _edge = 2
                    while not isinstance(self.parent_child[_node][2],EnumInfo):
                        if isinstance(self.parent_child[node][2],bool):
                            s.append(self.branch_transposition(_node,2,_node,0))
                        else:
                            _node = self.parent_child[_node][2]
                branches = self.branches()
                
                node1,edge1,dist1 = get_pair(self,node, edge, _node, _edge, branches)
#                 node2,edge2,dist2 = get_pair(self,_node, _edge, node, edge, branches) 
                if node1 != _node or edge1 != _edge:
                    
#                     if dist1 < dist2:
                    s.append(self.branch_transposition(_node,_edge,node1,edge1))
#                     else:
#                         s.append(self.branch_transposition(node,edge,node2,edge2))
            if edge == 2:
                n = node
                _node = self.parent_child[node].parent
                _edge = self.parent_child[_node].childs.index(node)
                while _edge == 2:
                    # Подъем до не Z edge
                    n = _node
                    _node = self.parent_child[_node].parent
                    _edge = sefl.parent_child[_node].childs.index(n)
                n = _node
                # cпуск по Z
                if _edge == 0:
                    if isinstance(self.parent_child[n][1],bool):
                        s.append(self.branch_transposition(n,1,n,2))
                        
                    if isinstance(self.parent_child[n][1],EnumInfo):
                        _node = n
                        _edge = 1

                    else:
                        _node = self.parent_child[n][1]
                        _edge = 2
                        while not isinstance(self.parent_child[_node][2],EnumInfo):
                            if isinstance(self.parent_child[_node][2],bool):
                                s.append(self.branch_transposition(_node,2,_node,1))
                            else:
                                _node = self.parent_child[_node][2]
                    
                    branches = self.branches()
                    node1,edge1,dist1 = get_pair(self,node, edge, _node, _edge, branches)
                    node2,edge2,dist2 = get_pair(self,_node, _edge, node, edge, branches) 
                    if node1 != _node or edge1 != _edge:

                        if dist1 < dist2:
                            s.append(self.branch_transposition(_node,_edge,node1,edge1))
                        else:
                            s.append(self.branch_transposition(node,edge,node2,edge2))
                if edge == 1:
                    if isinstance(self.parent_child[n][0],bool):
                        s.append(self.branch_transposition(n,0,n,2))
                    if isinstance(self.parent_child[n][0],EnumInfo):
                        _node = node
                        _edge = 0

                    else:
                        _node = self.parent_child[n][0]
                        _edge = 2
                        while not isinstance(self.parent_child[_node][2],EnumInfo):
                            if isinstance(self.parent_child[_node][2],bool):
                                s.append(self.branch_transposition(_node,2,_node,0))
                            else:
                                _node = self.parent_child[_node][2]
                    branches = self.branches()
                    node1,edge1,dist1 = get_pair(self,node, edge, _node, _edge, branches)
                    node2,edge2,dist2 = get_pair(self,_node, _edge, node, edge, branches) 
                    if node1 != _node or edge1 != _edge:

                        if dist1 < dist2:
                            s.append(self.branch_transposition(_node,_edge,node1,edge1))
                        else:
                            s.append(self.branch_transposition(node,edge,node2,edge2))
        return s
    
    def tojw(self):
        """
        Use less effective algorithms with transformation to JW tree
        """
        s = []
        def _find_first_branch(node):
            if self.parent_child[node][0]:
                return [node,0]
            if self.parent_child[node][1]:
                return [node,1]
            if self.parent_child[node][2]:
                return _find_first_branch(self.parent_child[node][2])
            return [False, False]
    # ищет первую не занятую z вершину   
        def _find_second_branch(node,edge):
            parents,flag = -1, False
            def down(parent):
                nonlocal parents, flag
                for index, child in enumerate(self.parent_child[parent].childs):
                    if flag:
                        break
                    if parent == node and index == edge:
                        i=0
                    elif isinstance(child,EnumInfo) and index == 2:                        
                        parents = parent
                        flag = True
                    elif child:
                        down(child)
            down(0)
            return parents
        def _check_xy():
            parentf = 0 
            indexf = 0
            def down(parent):
                nonlocal parentf,indexf 
                for index, child in enumerate(self.parent_child[parent].childs):
                    if isinstance(child, bool):
                        parentf, indexf = parent, index 
                    if child:
                        down(child)
            down(0)
            return [parentf,indexf]
        def check_pair(pair):
            if min(pair[0],pair[1]) % 2 == 1 and abs(pair[0] - pair[1]) == 1:
                return True
            else:
                return False
        def find_trans(pair_list,index):
            pair = pair_list[index]
            indexes = []
            if check_pair(pair):
                return False
            for i, ipair in enumerate(pair_list[index + 1:]):
                if check_pair([ipair[0],pair[0]]):
                    return [2*index + 1, 2*(i + index + 1) ]
                if check_pair([ipair[0],pair[1]]):
                    return [2*index, 2*(i + index + 1) ]
                if check_pair([ipair[1],pair[0]]):
                    return [2*index + 1, 2*(i + index + 1) + 1]
                if check_pair([ipair[1],pair[1]]):
                    return [2*index, 2*(i + index + 1) + 1]
                
        first, index = _check_xy()
#       Здесь сводится все к дереву JW
        s.append(self.branch_transposition(first, index, first, 2))
        s.append(self.branch_transposition(first, 0, _find_second_branch(first,0), 2))
        s.append(self.branch_transposition(first, 1, _find_second_branch(first,1), 2))
        while True:
            fl = _find_first_branch(0)
            if fl[0] or not isinstance(fl[0], bool):
                s.append(self.branch_transposition(fl[0], fl[1], _find_second_branch(fl[0],fl[1]), 2))
            else:
                break
#       Далее просто переставляются оставшиеся ветви
        branches, num_list = self.branches(True)
        num_list = [i.num for i in num_list]
        pair_list = [[num_list[2*i], num_list[2*i+1]] for i in range(len(num_list)//2)]
        L = len(pair_list)
        for index in range(L):
            k = find_trans(pair_list, index)
            if k:
                q1 = branches[k[0]][-1]
                q2 = branches[k[1]][-1]
                s.append(self.branch_transposition(q1[0], gate_index[q1[1]], q2[0], gate_index[q2[1]]))
            branches, num_list = self.branches(True)
            num_list = [i.num for i in num_list]
            pair_list = [[num_list[2*i], num_list[2*i+1]] for i in range(len(num_list)//2)]
        return s
    