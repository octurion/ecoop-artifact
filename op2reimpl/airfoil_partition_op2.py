#!/usr/bin/env python3
import sys
from itertools import chain, islice

def chunks(iterable, size=256):
    iterator = iter(iterable)
    for first in iterator:
        yield chain([first], islice(iterator, size - 1))

def edge_indices(f, num_edges):
    for i in range(0, num_edges):
        line = next(f).split()
        yield i, int(line[2]), int(line[3])

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

    graph_colors = []
    for chunk in chunks(edge_indices(f, num_edges)):
        min_idx = num_edges
        max_idx = 0
        indices = set()
        for i, e1, e2 in chunk:
            min_idx = min(min_idx, i)
            max_idx = max(max_idx, i)
            indices.add(e1)
            indices.add(e2)

        for r, s in graph_colors:
            if s.isdisjoint(indices):
                r.append((min_idx, max_idx + 1))
                s.update(indices)
                break
        else:
            new_color = [(min_idx, max_idx + 1)], indices
            graph_colors.append(new_color)

    for i, (r, _) in enumerate(graph_colors):
        print("Color: {}".format(i))
        for begin, end in r:
            print("({}, {})".format(begin, end))

    # sets = []
    # for e in edge_indices(f, num_edges):
    #     for s1, s2 in sets:
    #         if e[0] not in s1 and e[1] not in s2:
    #             s1.add(e[0])
    #             s2.add(e[1])
    #             break
    #     else:
    #         new_set = {e[0]}, {e[1]}
    #         sets.append(new_set)
    #
    # print([len(s1) for (s1, s2) in sets])

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        meme(f)
