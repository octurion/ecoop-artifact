#!/usr/bin/env python3
import csv
import itertools as it
import os
import sys

def main(root, machines):
    weight_variants = ["Aos", "Mixed", "Soa"]
    joint_variants = ["Scattered", "Pooled"]
    skip_first = 10
    num_runs = 10
    skip_every = 2

    for m in machines:
        cols = [
            ["1e5", "2e5", "3e5", "4e5", "5e5", "6e5", "7e5", "8e5", "9e5", "1e6"]
        ]
        filename = os.path.join(root, "stickmen_{}.csv".format(m))
        with open(filename, 'r') as f:
            r = csv.reader(f, delimiter=',')
            for _ in range(skip_first):
                next(r)

            for _ in range(len(weight_variants) * len(joint_variants)):
                new_col = []
                for _ in range(num_runs):
                    e = next(r)
                    cpu_time = e[3]
                    new_col.append(cpu_time)

                for _ in range(skip_every):
                    next(r)

                cols.append(new_col)

        filename = os.path.join(root, "stickmen_out_{}.csv".format(m))
        with open(filename, 'w') as f:
            w = csv.writer(f, delimiter=',')
            headers = ["{}-{}".format(j, i) for (i, j) in it.product(weight_variants, joint_variants)]
            headers.insert(0, "Count")
            w.writerow(headers)

            for e in zip(*cols):
                w.writerow(list(e))

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) >= 2 else os.getcwd()
    machines = sys.argv[2].split(",") if len(sys.argv) >= 3 else ["desktop", "laptop", "graphic", "ray", "voxel"]
    main(root, machines)
