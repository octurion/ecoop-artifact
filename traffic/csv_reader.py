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

def main(root):
    #machines = ["desktop", "laptop", "graphic", "ray", "voxel"]
    machines = ["laptop"]
    variants = ["aos", "mixed"]
    skip_lines = 10

    cols = []

    for m in machines:
        filename = os.path.join(root, "traffic_{}.csv".format(m))
        with open(filename, 'r') as f:
            new_cols = [[], []]

            r = csv.reader(f, delimiter=',')
            for _ in range(skip_lines):
                r.next()

            for e in r:
                test_name = e[0]
                cpu_time = e[3]

                idx = 0 if "Aos" in test_name else 1
                new_cols[idx].append(cpu_time)

            cols.extend(new_cols)

    filename = os.path.join(root, "traffic.csv")
    with open(filename, 'w') as f:
        write_csv(f, cols, machines, variants)

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) >= 2 else os.getcwd()
    main(root)
