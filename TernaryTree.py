from .BaseTree import BaseTree


gate_name = {0: 'X', 1: 'Y', 2: 'Z'}
gate_index = {'X': 0,'Y': 1, 'Z': 2}
dict_prod = {"II" : "I", "XX" : "I", "YY" : "I","ZZ" : "I", 
             "XY" : 'Z',"YX" : 'Z',"XZ" : 'Y',"ZX" : 'Y',"YZ" : 'X',"ZY" : 'X',
            "IX" : "X","XI" : "X","YI" : "Y","IY" : "Y","IZ" : "Z","ZI" : "Z"}

class TernaryTree(BaseTree):
    """
    Ternary tree object for initial state preparation
    """
    def to0vac(self):
        """
        Use algorithm described in the graduate work
        """
        help_dict = {0:1, 1:0}
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
       
        def get_num(node):
            while child:
                child = self.parent
                
        def get_Zpair(node,branches):
            child = node[0]
            parent = self.parent_child[node[0]].parent
            edge = self.parent_child[parent].childs.index(node[0]) 
            #подъем до не Z
            while edge == 2:
                child = parent
                parent = self.parent_child[parent].parent
                edge = self.parent_child[parent].childs.index(child) 
#             спуск до листа
            child = self.parent_child[parent][help_dict[edge]]
            while child:
                parent = child 
                child = self.parent_child[parent][2]
            return child.num, [parent, self.parent_child[parent].childs.index(child)]
        s = []
        
        first, index = _check_xy()
        parent = 0 
        child = self.parent_child[parent][2]
        while child:
            parent = child
            child = self.parent_child[parent][2]
        s.append(self.branch_transposition(first, index, parent,2))
#         s.append(self.branch_transposition(first, 0, _find_second_branch(first,0), 2))
#         s.append(self.branch_transposition(first, 1, _find_second_branch(first,1), 2))
        
        branches, num_list = self.branches(get_num = True)
#         num_list = [num.num for num in num_list]
        new_branches = [[] for _ in range(len(branches))]
        for i,branch in enumerate(branches):
            for j,node in enumerate(branch):
                (new_branches[i]).append([node[0],gate_index[node[1]]])
        branches = new_branches
        for index,branch in enumerate(branches):
            
            if branch[-1][1] == 0:
                if branches[index + 1][-1][1] == 1:
                    for _index,_node in enumerate([branch[-1] for branch in branches[index + 2:]]):
                        if check_pair([num_list[index],num_list[_index +index + 2]]):
                            s.append(self.branch_transposition(branches[index + 1][-1][0], 
                                    branches[index + 1][-1][1], branches[_index+ index + 2][-1][0], 
                                    branches[_index + index + 2][-1][1]))
                            num_list[index + 1],num_list[_index + index + 2] = num_list[_index + index + 2],num_list[index + 1]
                            break
            if branch[-1][1] == 2:
                num, _node = get_Zpair(branch[-1],branches)
                if not check_pair([num_list[index],num]):
                    for _index,_node in enumerate([branch[-1] for branch in branches]):
                        if check_pair([num,num_list[_index]]):
                            s.append(self.branch_transposition(branch[-1][0],branch[-1][1],branches[_index][-1][0], 
                                    branches[_index][-1][1]))
                            num_list[index ],num_list[_index] = num_list[_index],num_list[index ]
                            break
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
    