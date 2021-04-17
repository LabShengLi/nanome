FROM continuumio/miniconda3@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2
LABEL description="Nanocompare project" \
      author="yang.liu@jax.org"
RUN apt-get update -y \
  && apt-get install procps -y \
  && rm -rf /var/lib/apt/lists/*
  
# Create the environment:
COPY environment.yml /
RUN conda env create -f environment.yml && conda clean -a

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "nanocompare", "/bin/bash", "-c"]

ENV PATH /opt/conda/envs/cromwell-env/bin:$PATH
USER root
WORKDIR /data/
CMD ["bash"]
