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
    variants = ["one", "many"]
    skip_lines = 10

    prob50_cols = []
    prob70_cols = []
    prob90_cols = []

    for m in machines:
        filename = os.path.join(root, "doors_{}.csv".format(m))
        with open(filename, 'r') as f:
            prob50 = [[], []]
            prob70 = [[], []]
            prob90 = [[], []]

            r = csv.reader(f, delimiter=',')
            for _ in range(skip_lines):
                next(r)

            for e in r:
                test_name = e[0]
                cpu_time = e[3]

                idx = 0 if "OnePool" in test_name else 1

                if "probability:50" in test_name:
                    prob50[idx].append(cpu_time)
                elif "probability:70" in test_name:
                    prob70[idx].append(cpu_time)
                elif "probability:90" in test_name:
                    prob90[idx].append(cpu_time)
                else:
                    print("Warning: Couldn't parse line for test {}".format(test_name))

            prob50_cols.extend(prob50)
            prob70_cols.extend(prob70)
            prob90_cols.extend(prob90)

    filename = os.path.join(root, "doors50.csv")
    with open(filename, 'w') as f:
        write_csv(f, prob50_cols, machines, variants)

    filename = os.path.join(root, "doors70.csv")
    with open(filename, 'w') as f:
        write_csv(f, prob70_cols, machines, variants)

    filename = os.path.join(root, "doors90.csv")
    with open(filename, 'w') as f:
        write_csv(f, prob90_cols, machines, variants)

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) >= 2 else os.getcwd()
    machines = sys.argv[2].split(",") if len(sys.argv) >= 3 else ["desktop", "laptop", "graphic", "ray", "voxel"]
    main(root, machines)
