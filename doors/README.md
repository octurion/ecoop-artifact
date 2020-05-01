# Prerequisites

- Git
- CMake (version 3.5 or later)
- Any build system supported by CMake: Ninja (preferred and used in these
  instructions), Make, Visual Studio
- A compiler that supports C++11
- Any version of Python 3 to run the results script

# How to compile

- This repository depends on the Google Benchmarking library as a Git submodule.
  If you didn't populate the submodules when cloning the repository,
  run `git submodule update --init --recursive` in the root directory
- Run:
  1. `mkdir build && cd build`
  2. `cmake -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_TESTING=OFF -DBENCHMARK_ENABLE_GTEST_TESTS=OFF -DBENCHMARK_ENABLE_LTO=ON -GNinja ..`
  3. `ninja`

# How to run the benchmark

In the build directory (`build/`), run:

- `./bench --benchmark_out_format=csv --benchmark_out=doors_MACHINE_NAME.csv`
  (where `MACHINE_NAME` corresponds to the machine name, *i.e.,* `desktop`,
  `laptop`, `graphic`, `ray`, `voxel`).

# How to generate the result CSV files

- Consolidate all benchmark output CSV files in a single directory (*e.g.*
  `results/`)
- Then run `./csv_reader.py RESULTS_DIR` (where `RESULTS_DIR` is the
  aforementioned directory (*e.g.* `results/`).
- The output files should have the name `doors50.csv`, `doors70.csv`, and
  `doors90.csv`.

# Remarks

- If you want to run the CSV results generator with a different list of
  machine names, simply modify the `machines` list accordingly.
