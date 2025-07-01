# @Author   : Yang Liu
# @FileName : Dockerfile
# @Software : NANOME project
# @Organization : Li Lab
# @Website  : https://github.com/LabShengLi/nanome

# Set the base image to Ubuntu 18.04 and NVIDIA GPU from https://hub.docker.com/r/nvidia/cuda
# or from https://ngc.nvidia.com/catalog/containers/nvidia:cuda/tags
FROM nvidia/cuda:11.6.0-base-ubuntu18.04

# Author and maintainer
MAINTAINER Yang Liu <yang.liu@jax.org>
LABEL description="Nanome project in Sheng Li Lab" \
      author="liuyang2006@gmail.com"

# Guppy version 6.4.6 is not robust for meth call
# We choose v6.3.8 as stable one
ARG GUPPY_VERSION=6.3.8
ARG BUILD_PACKAGES="wget apt-transport-https procps git curl libnvidia-compute-460-server"
ARG DEBIAN_FRONTEND="noninteractive"

# Check version of guppy, ref: wget https://americas.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.1.3_linux64.tar.gz --no-check-certificate
# Another way: https://pypi.org/project/ont-pyguppy-client-lib/
# Install guppy-gpu version, ref: https://github.com/GenomicParisCentre/dockerfiles
RUN apt-get -q update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -q install --yes ${BUILD_PACKAGES} && \
    cd /tmp && \
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb  --no-check-certificate && \
    DEBIAN_FRONTEND="noninteractive" apt-get -q install --yes /tmp/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#Install miniconda
RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda.sh && \
    /bin/bash Miniconda.sh -b -p /opt/conda && \
    rm Miniconda.sh

# Adding conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Create the environment:
COPY guppy_stable_env.yml /
RUN conda env create --name guppy_stable --file=guppy_stable_env.yml && conda clean -a

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "guppy_stable", "/bin/bash", "-c"]

# Set env path into PATH
ENV PATH /opt/conda/envs/guppy_stable/bin:$PATH
ENV LD_LIBRARY_PATH /opt/conda/envs/guppy_stable/lib:$LD_LIBRARY_PATH
ENV HDF5_PLUGIN_PATH /opt/conda/envs/guppy_stable/hdf5/lib/plugin:$HDF5_PLUGIN_PATH

USER root
WORKDIR /data/

CMD ["bash"]
