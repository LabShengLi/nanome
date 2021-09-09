FROM genomicpariscentre/guppy-gpu:latest
LABEL description="Nanome project by Li Lab at JAX" \
      author="yang.liu@jax.org"
RUN apt-get update -y \
  && DEBIAN_FRONTEND=noninteractive apt-get install procps git curl -y \
  && rm -rf /var/lib/apt/lists/*

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common && \
    curl -O http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/cuda-ubuntu1604.pin && \
    mv cuda-ubuntu1604.pin /etc/apt/preferences.d/cuda-repository-pin-600 && \
    wget http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub && \
    apt-key add 7fa2af80.pub && \
    add-apt-repository "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/ /" && \
    apt update && \
    DEBIAN_FRONTEND=noninteractive apt -y install cuda

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
