import itertools as it

import csv
import os
import re
import statistics
import sys

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

def filter_op2_ours_file(f):
    for line in f:
        numlist = re.findall(r"^Wall clock time: (\d+.\d+)$", line)
        if numlist:
            yield numlist[0].strip()

def filter_op2_orig_file(f):
    for line in f:
        numlist = re.findall(r"^Total plan time:   (\d+.\d+)$", line)
        if numlist:
            yield numlist[0].strip()

        numlist = re.findall(r"^Max total runtime = (\d+.\d+)$", line)
        if numlist:
            yield numlist[0].strip()

def write_results(f, machines, variants, results):
    w = csv.writer(f, delimiter=',')
    headers = ["{}-{}".format(m, v) for m, v in it.product(machines, variants)]
    w.writerow(headers)
    for e in zip(*results):
        w.writerow(list(e))

def main(orig_root, ours_root, out_root):
    num_tries = 20

    #machines = ["desktop", "laptop", "graphic", "ray", "voxel"]
    machines = ["laptop"]
    airfoil_variants = ["orig", "aos", "mixed", "soa"]
    aero_variants = ["orig", "aos", "mixed"]

    aero_results = []
    airfoil_results = []

    for machine in machines:
        op2_orig = "op2_orig_{}.txt".format(machine)
        with open(os.path.join(orig_root, op2_orig), 'r') as f:
            lines = [e for e in filter_op2_orig_file(f)]
            plan_times = lines[0::2]
            execution_times = lines[1::2]

            times = [float(b) - float(a) for (a, b) in zip(plan_times, execution_times)]
            lists = [l for l in chunks(times, num_tries)]

            airfoil_results.append(lists[0])
            aero_results.append(lists[1])

        op2_ours = "op2_ours_{}.txt".format(machine)
        with open(os.path.join(ours_root, op2_ours), 'r') as f:
            lines = [e.strip() for e in filter_op2_ours_file(f)]
            lists = [l for l in chunks(lines, num_tries)]

            aero_results.extend(lists[:2])
            airfoil_results.extend(lists[2:])

    with open(os.path.join(out_root, 'airfoil.csv'), 'w') as f:
        write_results(f, machines, airfoil_variants, airfoil_results)

    with open(os.path.join(out_root, 'aero.csv'), 'w') as f:
        write_results(f, machines, aero_variants, aero_results)

#    for l in lists:
#        l = [float(e) for e in l]
#        print("Mean: {}".format(statistics.mean(l)))
#        print("Median: {}".format(statistics.median(l)))
#        print("Stdev: {}".format(statistics.stdev(l)))
#        print("Min: {}".format(min(l)))
#        print("Max: {}".format(max(l)))

if __name__ == "__main__":
    orig_path = sys.argv[1] if len(sys.argv) >= 2 else os.getcwd()
    ours_path = sys.argv[2] if len(sys.argv) >= 3 else os.getcwd()
    out_path = os.getcwd()
    main(orig_path, ours_path, out_path)

