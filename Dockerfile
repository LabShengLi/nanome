FROM genomicpariscentre/guppy-gpu:4.2.2
LABEL description="Nanocompare project" \
      author="yang.liu@jax.org"
RUN apt-get update -y \
  && apt-get install procps -y \
  && rm -rf /var/lib/apt/lists/*
 
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
SHELL ["conda", "run", "-n", "nanocompare", "/bin/bash", "-c"]

ENV PATH /opt/conda/envs/nanocompare/bin:$PATH
USER root
WORKDIR /data/
CMD ["bash"]
