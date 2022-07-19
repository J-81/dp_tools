FROM ubuntu:20.04

# Ensure no user interaction is requested
ARG DEBIAN_FRONTEND=noninteractive

# create group and user and install packages
RUN groupadd -r genuser && \
    useradd -g genuser genuser && \
    mkdir /home/genuser && \
    chmod -R 777 /home/genuser && \
    apt-get update && \
    apt-get install software-properties-common -y && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get install python3.10 python3.10-distutils curl samtools -y && \
    # ensure python3.10 is linked to python for shebang support
    ln -s /usr/bin/python3.10 /usr/bin/python
    

# copy dp_tools into container
COPY . /app

# set ownership and permissions
RUN chown -R genuser:genuser /app && \
    curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10 && \
    pip install /app &&  \
    # save space in image by removing source code after pip install
    rm -rf /app && \
    echo "export PATH=/home/genuser/.local/bin:$PATH" >> ~/.bashrc

# swith to user
USER genuser

WORKDIR /home/genuser