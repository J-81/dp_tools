FROM nfcore/gitpod:latest

# Ensure no user interaction is requested
ARG DEBIAN_FRONTEND=noninteractive

USER gitpod

# Install conda and python packages
RUN conda install -c conda-forge -c bioconda python=3.10 pytest samtools && \
    pip install mypy black
    
