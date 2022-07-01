FROM gitpod/workspace-full:latest

USER gitpod

RUN pyenv install miniconda3-4.7.12 && \
    pyenv activate miniconda3-4.7.12 && \
    conda env create -f condaEnv.yaml