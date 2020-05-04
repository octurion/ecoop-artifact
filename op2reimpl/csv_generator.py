#!/usr/bin/env python3
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

def main(root, machines):
    num_tries = 20

    airfoil_variants = ["orig", "aos", "mixed", "soa"]
    aero_variants = ["orig", "aos", "mixed"]

    aero_results = []
    airfoil_results = []

    for machine in machines:
        op2_orig = "op2_orig_{}.txt".format(machine)
        with open(os.path.join(root, op2_orig), 'r') as f:
            lines = [e for e in filter_op2_orig_file(f)]
            plan_times = lines[0::2]
            execution_times = lines[1::2]

            times = [float(b) - float(a) for (a, b) in zip(plan_times, execution_times)]
            lists = [l for l in chunks(times, num_tries)]

            airfoil_results.append(lists[0])
            aero_results.append(lists[1])

        op2_ours = "op2_ours_{}.txt".format(machine)
        with open(os.path.join(root, op2_ours), 'r') as f:
            lines = [e.strip() for e in filter_op2_ours_file(f)]
            lists = [l for l in chunks(lines, num_tries)]

            aero_results.extend(lists[:2])
            airfoil_results.extend(lists[2:])

    with open(os.path.join(root, 'airfoil.csv'), 'w') as f:
        write_results(f, machines, airfoil_variants, airfoil_results)

    with open(os.path.join(root, 'aero.csv'), 'w') as f:
        write_results(f, machines, aero_variants, aero_results)

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) >= 2 else os.getcwd()
    machines = sys.argv[2].split(",") if len(sys.argv) >= 3 else ["desktop", "laptop", "graphic", "ray", "voxel"]
    main(root, machines)
