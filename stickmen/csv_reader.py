import csv
import itertools as it
import os
import sys

def main(root):
    #machines = ["desktop", "laptop", "graphic", "ray", "voxel"]
    machines = ["laptop"]
    weight_variants = ["Aos", "Mixed", "Soa"]
    joint_variants = ["Scattered", "Pooled"]
    skip_first = 9
    num_runs = 10
    skip_every = 2

    for m in machines:
        cols = []
        filename = os.path.join(root, "stickmen_{}.csv".format(m))
        with open(filename, 'r') as f:
            r = csv.reader(f, delimiter=',')
            for _ in range(skip_first):
                r.next()

            for _ in range(len(weight_variants) * len(joint_variants)):
                new_col = []
                for _ in range(num_runs):
                    e = r.next()
                    cpu_time = e[3]
                    new_col.append(cpu_time)

                for _ in range(skip_every):
                    r.next()

                cols.append(new_col)

        filename = os.path.join(root, "stickmen_out_{}.csv".format(m))
        with open(filename, 'w') as f:
            w = csv.writer(f, delimiter=',')
            headers = ["{}-{}".format(j, i) for (i, j) in it.product(weight_variants, joint_variants)]
            w.writerow(headers)

            for e in zip(*cols):
                w.writerow(list(e))

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) >= 2 else os.getcwd()
    main(root)
