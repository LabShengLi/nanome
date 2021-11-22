# @Author   : Yang Liu
# @FileName : Dockerfile
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

# Set the base image to Ubuntu 18.04 and NVIDIA GPU from https://hub.docker.com/r/nvidia/cuda
# or from https://ngc.nvidia.com/catalog/containers/nvidia:cuda/tags
FROM nvidia/cuda:11.4.1-base-ubuntu18.04

# Author and maintainer
MAINTAINER Yang Liu <yang.liu@jax.org>
LABEL description="Nanome project in Li Lab at The Jackson Laboratory" \
      author="yang.liu@jax.org"

# Guppy version
ARG GUPPY_VERSION=5.0.16
ARG BUILD_PACKAGES="wget apt-transport-https procps git curl"
ARG DEBIAN_FRONTEND="noninteractive"
ARG NANOME_DIR=/opt/nanome
ARG METEORE_GITHUB="https://github.com/comprna/METEORE/archive/refs/tags/v1.0.0.tar.gz"
ARG DEEPSIGNAL_MODEL="https://zenodo.org/record/5513090/files/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7%2B.tar.gz"
ARG MEGALODON_MODEL="https://zenodo.org/record/5513090/files/megalodon_model.tar.gz"
ARG ECOLI_GENOME="https://zenodo.org/record/5513090/files/ecoli.tar.gz"
ARG HG38_CHR22_GENOME="https://zenodo.org/record/5513090/files/hg38_chr22.tar.gz"

ARG MEGALODON_VERSION=2.3.4

# Install guppy-gpu version, ref: https://github.com/GenomicParisCentre/dockerfiles
RUN apt-get -q update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -q install --yes ${BUILD_PACKAGES} libnvidia-compute-460-server && \
    cd /tmp && \
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${GUPPY_VERSION}-1~bionic_amd64.deb && \
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
COPY environment.yml /
RUN conda env create --name nanome --file=environment.yml && conda clean -a

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "nanome", "/bin/bash", "-c"]

# Install latest version for megalodon, even conflicts with fast5mod, they can work
RUN pip install megalodon==${MEGALODON_VERSION} && pip cache purge
RUN npm install -g inliner && npm cache clean --force

# Set nanome env path into PATH
ENV PATH /opt/conda/envs/nanome/bin:$PATH
USER root
WORKDIR /data/

# Get METEORE dir into /data/METEORE-1.0.0,
# Get DeepSignal model into /data/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+
# Get E. coli and hg38 chr22 regenome into /data
RUN cd /data && wget -q ${METEORE_GITHUB} &&\
    tar -xzf v1.0.0.tar.gz &&\
    rm -f v1.0.0.tar.gz &&\
    rm -rf METEORE-1.0.0/example_results METEORE-1.0.0/data &&\ 
    wget -q ${DEEPSIGNAL_MODEL} &&\
    wget -q ${MEGALODON_MODEL} &&\
    wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig &&\
    chmod +x bedGraphToBigWig &&\
    mv bedGraphToBigWig  /usr/local/bin/

CMD ["bash"]
