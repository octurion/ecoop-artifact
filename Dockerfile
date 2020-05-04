FROM ubuntu:16.04

COPY . /root
RUN apt-get update
RUN apt-get -y install --no-install-recommends git python3 build-essential cmake g++ libgomp1 ninja-build texlive-fonts-recommended texlive-pictures xz-utils
