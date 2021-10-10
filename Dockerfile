# Set the base image to Ubuntu 18.04 and NVIDIA GPU from https://hub.docker.com/r/nvidia/cuda
# or from https://ngc.nvidia.com/catalog/containers/nvidia:cuda/tags
FROM nvidia/cuda:11.4.1-base-ubuntu18.04

# Author and maintainer
MAINTAINER Yang Liu <yang.liu@jax.org>
LABEL description="Nanome project in Li Lab at The Jackson Laboratory" \
      author="yang.liu@jax.org"

# Guppy version
ARG PACKAGE_VERSION=5.0.14
ARG BUILD_PACKAGES="wget apt-transport-https"
ARG DEBIAN_FRONTEND=noninteractive
ARG NANOME_DIR=/opt/nanome

# Install guppy-gpu version, ref: https://github.com/GenomicParisCentre/dockerfiles
RUN apt-get -q update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -q install --yes $BUILD_PACKAGES libnvidia-compute-460-server && \ 
    cd /tmp && \
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${PACKAGE_VERSION}-1~bionic_amd64.deb && \
    DEBIAN_FRONTEND="noninteractive" apt-get -q install --yes /tmp/ont_guppy_${PACKAGE_VERSION}-1~bionic_amd64.deb && \
    rm *.deb && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install useful tools, e.g., git, curl, etc.
RUN apt-get -q update -y &&\
    DEBIAN_FRONTEND="noninteractive" apt-get install procps git curl -y &&\
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
COPY environment.yml /
RUN conda env create --name nanome --file=environment.yml && conda clean -a

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "nanome", "/bin/bash", "-c"]

# Install latest version for megalodon, even conflicts with fast5mod, they can work
RUN pip install megalodon

# Set nanome env path into PATH
ENV PATH /opt/conda/envs/nanome/bin:$PATH
USER root
WORKDIR /data/


RUN mkdir -p ${NANOME_DIR}

# Copy additonal scripts
ADD inputs/ ${NANOME_DIR}/inputs
ADD Rshiny/ ${NANOME_DIR}/Rshiny
ADD src/ ${NANOME_DIR}/src
ADD test_data/ ${NANOME_DIR}/test_data
ADD utils/ ${NANOME_DIR}/utils

# Copy nextflow scripts
ADD conf/ ${NANOME_DIR}/conf
ADD main.nf ${NANOME_DIR}/
ADD nextflow.config ${NANOME_DIR}/
ADD README.md ${NANOME_DIR}/
ADD LICENSE ${NANOME_DIR}/

# Change execute permissions
RUN find ${NANOME_DIR} -name "*.py" -type f -exec chmod +x {} \;
RUN find ${NANOME_DIR} -name "*.sh" -type f -exec chmod +x {} \;

CMD ["bash"]
