from . import EnumInfo, NodeInfo
def get_bin_parent_child(n):
    if n == 4:
        parent_child = {0 : NodeInfo(-1,[1,2,False])}
        parent_child[1] = NodeInfo(0,[3,EnumInfo(4), False])
        parent_child[2] = NodeInfo(0,[4, EnumInfo(8), False])
        parent_child[3] =  NodeInfo(1,[EnumInfo(1),EnumInfo(2),EnumInfo(3)])
        parent_child[4] = NodeInfo(2,[EnumInfo(5),EnumInfo(6),EnumInfo(7)])
        return parent_child
    if n == 6:
        parent_child = {0 : NodeInfo(-1,[1,2,False])}
        parent_child[1] = NodeInfo(0,[3,4, False])
        parent_child[2] = NodeInfo(0,[5,6, False])
        parent_child[3] =  NodeInfo(1,[EnumInfo(1),EnumInfo(2),EnumInfo(3)])
        parent_child[4] = NodeInfo(1,[EnumInfo(4),EnumInfo(5),EnumInfo(6)])
        parent_child[5] =  NodeInfo(2,[EnumInfo(7),EnumInfo(8),EnumInfo(9)])
        parent_child[6] =  NodeInfo(2,[EnumInfo(10),EnumInfo(11),EnumInfo(12)])
        return parent_child
    if n == 8:
        parent_child = {0 : NodeInfo(-1,[1,2,False])}
        parent_child[1] = NodeInfo(0,[3,4,5])
        parent_child[2] = NodeInfo(0,[6,7, 8])
        parent_child[3] = NodeInfo(1,[EnumInfo(1),EnumInfo(2), EnumInfo(3)])
        parent_child[4] = NodeInfo(1,[EnumInfo(4),EnumInfo(5), EnumInfo(6)])
        parent_child[5] = NodeInfo(1,[EnumInfo(7),EnumInfo(8), False])
        parent_child[6] = NodeInfo(2,[EnumInfo(9),EnumInfo(10), EnumInfo(11)])
        parent_child[7] = NodeInfo(2,[EnumInfo(12),EnumInfo(13),EnumInfo(14)])
        parent_child[8] = NodeInfo(2,[EnumInfo(15),EnumInfo(16),False])

        return parent_child
#     if n == 6:
#         return parent_child
#     if n == 6:
#         return parent_child