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
the source code using Docker:

1. Run `docker build .` the following in the top level directory
(i.e. where the `Dockerfile` resides):

- `docker build .`

Approximately 150 MB of APT dependencies will be downloaded at first.
After APT finishes up, building should take approximately 2-3 minutes.

2. When building finishes, you will see a message along the terms of:

`Successfully built cfe3ef41e886`

Make note of this hexadecimal container identifier (let's call it
`YOUR_CONTAINER_ID`), you'll need it to run the case studies.

## Without Docker

1. We assume you are running Ubuntu 16.04. We have not tested if there are any
   additional dependencies that need to be installed on newer Ubuntu versions.

2. Install the necessary dependencies by running:

- `sudo apt-get -y install --no-install-recommends git cmake g++ libgomp1 ninja-build texlive-fonts-recommended texlive-pictures xz-utils`

3. Run `./build.sh`

# How to run the case studies

In our paper, we run the case studies on multiple machines; we then consolidate
the results into one directory and process them to generate unified results
files.

We cannot expect you to do the same. As such, we also provide a method to
generate the results on just a single machine.

Notice that running the case studies will take approximately 12 hours, so
make sure to let the case studies run overnight :)

## If you want to run the case studies on a single machine

1. If you are using Docker, run:

- `mkdir local_mount`
- `docker run --mount type=bind,source=${PWD}/local_mount,destination=/mount -w /root -it YOUR_CONTAINER_ID`

Where `YOUR_CONTAINER_ID` is the container identifier from above.

While inside the Docker container, you can transfer files from the Docker
container's filesystem to the local filesystem by copying them to the `/mount`
container directory. Any copied files will be accessible from the local file
system in the `./local_mount/` directory.

2. Run the following:

- `./measure.sh -o=results/ -m=single all`
- `./gather.sh -m=single -i=results/ -o=charts/csv_data/ all`
- `./gen_charts.sh --single`

3. The resulting charts are located in `build/charts/charts_single.pdf`

4. To extract the PDF from the Docker container:

- **Don't exit the container yet!**
- Run `cp build/charts/charts_single.pdf /mount`

The chart should be located in `./local_mount/charts_single.pdf`

## If you want to run the case studies on multiple machines

1. Come up with names for your machines. We use `desktop`, `laptop`, `graphic`,
`ray`, `voxel` as our machine names in our paper.

2. If you are using Docker, on each machine, run:

- `mkdir local_mount`
- `docker run --mount type=bind,source=${PWD}/local_mount,destination=/mount -w /root -it YOUR_CONTAINER_ID`

Where `YOUR_CONTAINER_ID` is the container identifier from above.

While inside the Docker container, you can transfer files from the Docker
container's filesystem to the local filesystem by copying them to the `/mount`
container directory. Any copied files will be accessible from the local file
system in the `./local_mount/` directory.

3. On each machine, run:

- `./measure.sh -o=results/ -m=_MACHINE_NAME_HERE_ all`

4. Manual labour: Consolidate all generated `csv` and `txt` files from the
  `results/` directories on all machines in one directory (e.g.
  `consolidated/`).

  If using Docker, you can copy the generated files to the `/mount` directory.
  The generated files will be accessible from `./local_mount/`

5. Run the following:

- `./gather.sh -m=desktop,laptop,graphic,ray,voxel -i=consolidated/ -o=charts/csv_data/ all`
- `./gen_charts.sh`

6. The resulting charts are located in `build/charts/charts.pdf`

7. To extract the PDF from the Docker container:

- **Don't exit the container yet!**
- Run `cp build/charts/charts.pdf /mount`

The chart should be located in `./local_mount/charts.pdf`
