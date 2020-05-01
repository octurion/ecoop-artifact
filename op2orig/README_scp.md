# Prerequisites

- Git
- CMake (version 3.5 or later)
- Any build system supported by CMake: Ninja (preferred and used in these
  instructions), Make, Visual Studio
- A compiler that supports C++11 and OpenMP 2.0 or later (*e.g.* `gcc` via
  the `libgomp` library)
- GNU Octave
- Any version of Python 3 to run the results script

# Modifications

Our modifications on the original code consist of the following:
- We eliminated the use of a nonstandard `gcc` C++ extension (empty structs),
  which was causing us compiler errors in our version of `gcc`. The
  `p0.patch` file consists of the relevant changes (these can also be seen
  through `git diff HEAD~1`.
- We added a `script.sh` file in the root directory for automating the
  generation of the results.

# How to compile

1. Run:
   - `cd apps/c/`
   - `mkdir build/`
   - `cp ../../script.sh build/`
2. The two OP2 applications rely on two input files named `new_grid.dat` and
  `FE_grid.dat` (generated inside the `build` directory). These are generated
  via GNU Octave and they will take a significant amount of time to
  be generated. If you do not wish to wait for them to be generated, you can
  download them from here:
  - [new_grid.dat](https://drive.google.com/open?id=1Afw2IDrObOQK7-T-pjCX30RWj-1gL3et)
  - [FE_grid.dat](https://drive.google.com/open?id=144ydk9w5oehzERWHq68SmlBtOpcecr_0)
  And place `new_grid.dat` and `FE_grid.dat` inside the `build/` directory.
3. Run `./cmake.local -DCMAKE_BUILD_TYPE=Release`

# How to run the benchmark

In the root directory, run:

- `./script.sh > op2_orig_MACHINE_NAME.txt`
  (where `MACHINE_NAME` corresponds to the machine name, *i.e.,* `desktop`,
  `laptop`, `graphic`, `ray`, `voxel`).

# How to generate the result CSV files

- Consolidate all benchmark output `txt` files corresponding to
  the original OP2 applications' runs (*i.e.* repo `OP2-Common`)
  on the various machines, as well as the output `csv` files
  corresponding to the reimplementation's runs (*i.e.* repo
  `op2reimpl`) in a single directory (*e.g.* `results/`).
- Run `./csv_generator.py RESULTS_DIR` (where `RESULTS_DIR` is the
  aforementioned directory (*e.g.* `results/` and `csv_generator.py`
  is located in repo `op2reimpl`).
- The output files should have the name `aero.csv` and `airfoil.csv`.
