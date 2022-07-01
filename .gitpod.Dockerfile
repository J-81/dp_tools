FROM continuumio/miniconda3

# Install conda environment for this project
RUN wget https://github.com/J-81/dp_tools/raw/gitpod_setups/condaEnv.yaml && \
    cat condaEnv.yaml && \
    conda env create -f condaEnv.yaml