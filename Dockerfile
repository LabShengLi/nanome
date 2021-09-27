# Set the base image to Ubuntu 18.04 and NVIDIA GPU from https://hub.docker.com/r/nvidia/cuda
FROM nvidia/cuda:11.4.1-base-ubuntu18.04

# Author and maintainer
MAINTAINER Yang Liu <yang.liu@jax.org>
LABEL description="Nanome project in Li Lab at The Jackson Laboratory" \
      author="yang.liu@jax.org"

# Guppy version
ARG PACKAGE_VERSION=5.0.14
ARG BUILD_PACKAGES="wget apt-transport-https"
ARG DEBIAN_FRONTEND=noninteractive

# Install guppy-gpu version, ref: https://github.com/GenomicParisCentre/dockerfiles
RUN apt update && \
    apt install --yes $BUILD_PACKAGES libnvidia-compute-460-server && \ 
    cd /tmp && \
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${PACKAGE_VERSION}-1~bionic_amd64.deb && \
    apt install --yes /tmp/ont_guppy_${PACKAGE_VERSION}-1~bionic_amd64.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install useful tools, e.g., git, curl, etc.
RUN apt-get update -y \
  && DEBIAN_FRONTEND=noninteractive apt-get install procps git curl -y \
  && rm -rf /var/lib/apt/lists/*

# Install cuda
# RUN apt-get update -y && \
#     DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common && \
#     curl -O http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/cuda-ubuntu1604.pin && \
#     mv cuda-ubuntu1604.pin /etc/apt/preferences.d/cuda-repository-pin-600 && \
#     wget http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub && \
#     apt-key add 7fa2af80.pub && \
#     add-apt-repository "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/ /" && \
#     apt update && \
#     DEBIAN_FRONTEND=noninteractive apt -y install cuda

#Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda.sh && \
    /bin/bash Miniconda.sh -b -p /opt/conda && \
    rm Miniconda.sh

# Adding conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Create the environment:
COPY environment.yml /
RUN conda env create -f environment.yml && conda clean -a

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "nanome", "/bin/bash", "-c"]

# Install latest version for megalodon, even conflicts with fast5mod, they can work
RUN pip install megalodon

# Set nanome env path into PATH
ENV PATH /opt/conda/envs/nanome/bin:$PATH
USER root
WORKDIR /data/

# Copy additonal scripts
RUN mkdir /opt/bin
ADD utils/ /opt/bin/utils
ADD src/ /opt/bin/src
ADD test_data/ /opt/bin/test_data
ADD inputs/ /opt/bin/inputs
RUN chmod +x /opt/bin/utils/*
RUN chmod +x /opt/bin/src/*
RUN chmod +x /opt/bin/src/nanocompare/*
ENV PATH="$PATH:/opt/bin/utils"
ENV PATH="$PATH:/opt/bin/src"
ENV PATH="$PATH:/opt/bin/src/nanocompare"
ENV PYTHONPATH="/opt/bin/src"

CMD ["bash"]
