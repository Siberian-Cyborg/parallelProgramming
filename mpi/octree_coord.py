def node_idx(c):
    if len(c) == 0:
        return 0

    idx = 0
    k = len(c)
    for i in range(k):
        idx += pow(8, i)

    assert idx == ((1<<(k*3))-1)/7

    for i,v in enumerate(reversed(c)):
        idx += pow(8,i) * v

    return idx

if False:
    print (node_idx([]))
    print (node_idx([0]))
    print (node_idx([7]))
    print (node_idx([0,0,0]))
    print (node_idx([0,0,1]))
    print (node_idx([0,0,7]))
    print (node_idx([0,1,0]))
    print (node_idx([0,7,7]))
    print (node_idx([1,0,0]))
    print (node_idx([7,7,7]))
    print (node_idx([0,0,0,0]))

    for i in range(8):
        print(node_idx([i]))

    for i in range(8):
        for j in range(8):
            print(node_idx([i, j]))

    for i in range(8):
        for j in range(8):
            for k in range(8):
                print(node_idx([i, j, k]))


def construct_tree(max_depth, cur_depth, c):
    print (c, node_idx(c))
    if cur_depth < max_depth:
        for i in range(8):
            construct_tree(max_depth, cur_depth+1, c + [i])

construct_tree(3, 0, [])
