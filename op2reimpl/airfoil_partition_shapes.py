#!/usr/bin/env python3
import sys

def edge_indices(f, num_edges):
    for i in range(0, num_edges):
        line = next(f).split()
        yield line[2], line[3]

def bedge_indices(f, num_edges):
    for i in range(0, num_edges):
        line = next(f).split()
        yield line[2]

def meme(f):
    header = next(f).split()
    num_nodes = int(header[0])
    num_cells = int(header[1])
    num_edges = int(header[2])
    num_bedges = int(header[3])

    print("Num of nodes: {}".format(num_nodes))
    print("Num of cells: {}".format(num_cells))
    print("Num of edges: {}".format(num_edges))
    print("Num of backedges: {}".format(num_bedges))

    for _ in range(0, num_nodes + num_cells):
        next(f)

    sets = []
    for e in edge_indices(f, num_edges):
        for s1, s2 in sets:
            if e[0] not in s1 and e[1] not in s2:
                s1.add(e[0])
                s2.add(e[1])
                break
        else:
            new_set = {e[0]}, {e[1]}
            sets.append(new_set)

    print([len(s1) for (s1, s2) in sets])

    dank_sets = []
    for e in bedge_indices(f, num_edges):
        for s in dank_sets:
            if e not in s:
                s.add(e)
                break
        else:
            new_set = {e}
            dank_sets.append(new_set)

    print([len(s) for s in dank_sets])

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        meme(f)
