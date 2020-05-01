#!/usr/bin/env python3
import sys

def cell_indices(f, num_cells):
    for i in range(0, num_cells):
        line = [int(e) for e in next(f).split()]
        yield line

def meme(f):
    header = next(f).split()
    num_nodes = int(header[0])
    num_cells = int(header[1])
    num_bnodes = int(header[2])

    print("Num of nodes: {}".format(num_nodes))
    print("Num of cells: {}".format(num_cells))
    print("Num of backnodes: {}".format(num_bnodes))

    for _ in range(0, num_nodes):
        next(f)

    sets = []
    for n, e in enumerate(cell_indices(f, num_cells)):
        for v, s in sets:
            if all(i not in s for i in e):
                v.append(n)
                s.update(e)
                break
        else:
            sets.append(([n], set(e)))

    print([len(v) for v, s in sets])

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        meme(f)
