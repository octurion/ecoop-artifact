FROM ubuntu:16.04

RUN apt-get update \
    && apt-get -y install --no-install-recommends \
        git \
        python3 \
        build-essential \
        cmake \
        g++ \
        libgomp1 \
        ninja-build \
        texlive-fonts-recommended \
        texlive-pictures \
        xz-utils \
    && rm -rf /var/lib/apt/lists/*

COPY . /root

WORKDIR /root
RUN ./build.sh
