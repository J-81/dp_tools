FROM ubuntu:20.04

# Ensure no user interaction is requested
ARG DEBIAN_FRONTEND=noninteractive

# create group and user and install packages
RUN groupadd -r genuser && \
    useradd -g genuser genuser && \
    mkdir /home/genuser && \
    chmod -R 777 /home/genuser && \
    apt-get update
RUN apt-get install software-properties-common -y && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get install python3.10 python3.10-distutils curl -y
    

# copy dp_tools into container
COPY . /app

# set ownership and permissions
RUN chown -R genuser:genuser /app

# swith to user
USER genuser

RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10 && \
    /home/genuser/.local/bin/pip install /app &&  \
    echo "export PATH=/home/genuser/.local/bin:$PATH" >> ~/.bashrc

WORKDIR /home/genuser