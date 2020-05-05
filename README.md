# Introduction

This is the code corresponding to the case studies in our paper. It is
distributed into two variants: As a Git repository and as the artifact
submitted to our paper. This README provides instructions on how to build
the source code in both variants, run the case studies, collect the results,
and generate the charts in our paper.

# Before you start

Our code depends on the
[Google C++ Benchmark library](https://github.com/google/benchmark).
If you downloaded the source code as an artifact, then you can simply skip to
*How to build* below; the source code of the library is already included in the
artifact.

If you are cloning the repository, make sure to run `git clone` with the
`--recurse-submodules` flag (or run `git submodule update --init --recursive`
afterwards).

# How to build

## With Docker

We recommend that you use [Docker](https://docker.io) for building. To build
the source code using Docker, run the following in the top level directory
(i.e. where the `Dockerfile` resides):

- `docker build .`
- `docker run -it CONTAINER_HASH` (where `CONTAINER_HASH` is the hash of the
  newly generated Docker container).
- `cd /root`
- `./build.sh`

Approximately 150 MB of APT dependencies will be downloaded at first.
After APT finishes up, building should take approximately 2-3 minutes.

## Without Docker

On Ubuntu 16.04, you will need to install the necessary dependencies by running:

`sudo apt-get -y install --no-install-recommends git cmake g++ libgomp1 ninja-build texlive-fonts-recommended texlive-pictures xz-utils`

We have not tested if there are any additional dependencies that need to be
installed on newer Ubuntu versions.

Afterwards, run `./build.sh`

# How to run the case studies

In our paper, we run the case studies on multiple machines; we then consolidate
the results into one directory and process them to generate unified results
files.

We cannot expect you to do the same. As such, we also provide a method to
generate the results on just a single machine.

Notice that running the case studies will take approximately 12 hours, so
make sure to let the case studies run overnight :)

## If you want to run the case studies on a single machine

In the root directory (i.e. the directory of the cloned Git repository or
the Docker `/root/` directory), run the following:

- `./measure.sh -o=results/ -m=single all`
- `./gather.sh -m=single -i=results/ -o=charts/csv_data/ all`
- `./gen_charts.sh --single`
- The resulting charts are located in `build/charts/charts_single.pdf`

## If you want to run the case studies on multiple machines

Come up with names for your machines. We use `desktop`, `laptop`, `graphic`,
`ray`, `voxel` as our machine names.

- On each machine, run: `./measure.sh -o=results/ -m=_MACHINE_NAME_HERE_ all`
- Manual labour: Consolidate all generated `csv` and `txt` files from the
  `results/` directories on all machines in one directory (e.g.
  `consolidated/`).
- `./gather.sh -m=single -i=consolidated/ -o=charts/csv_data/ all`
- `./gen_charts.sh`
- The resulting charts are located in `build/charts/charts.pdf`
