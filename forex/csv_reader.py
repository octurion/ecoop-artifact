#!/usr/bin/env python3
import csv
import itertools as it
import os
import sys

def write_csv(f, cols, machines, variants):
    w = csv.writer(f, delimiter=',')
    headers = ["{}-{}".format(m, v) for (m, v) in it.product(machines, variants)]
    w.writerow(headers)

    for e in zip(*cols):
        w.writerow(list(e))

def main(root, machines):
    variants = ["one-pool", "many-pools", "many-pools-soa"]
    skip_lines = 10

    cols = []

    for m in machines:
        filename = os.path.join(root, "forex_{}.csv".format(m))
        with open(filename, 'r') as f:
            new_cols = [[], [], []]

            r = csv.reader(f, delimiter=',')
            for _ in range(skip_lines):
                next(r)

            for e in r:
                test_name = e[0]
                cpu_time = e[3]

                if "OnePool/" in test_name:
                    idx = 0
                elif "ManyPools/" in test_name:
                    idx = 1
                elif "ManyPoolsSoa/" in test_name:
                    idx = 2
                else:
                    continue

                new_cols[idx].append(cpu_time)

            cols.extend(new_cols)

    filename = os.path.join(root, "forex.csv")
    with open(filename, 'w') as f:
        write_csv(f, cols, machines, variants)

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) >= 2 else os.getcwd()
    machines = sys.argv[2].split(",") if len(sys.argv) >= 3 else ["desktop", "laptop", "graphic", "ray", "voxel"]
    main(root, machines)
