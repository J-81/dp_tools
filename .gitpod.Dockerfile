FROM gitpod/workspace-full:2022-07-26-05-11-39

# Ensure no user interaction is requested
ARG DEBIAN_FRONTEND=noninteractive

USER gitpod

# create group and user and install packages
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b && \
    export PATH=/home/gitpod/miniconda3/bin:$PATH && \
    echo "PATH=/home/gitpod/miniconda3/bin:$PATH" >> /home/gitpod/.bashrc && \
    conda install -c conda-forge -c bioconda python=3.10 pytest samtools && \
    pip install mypy black
    
