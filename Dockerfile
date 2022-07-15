FROM python:3.10.5

# create group and user
RUN groupadd -r genuser && \
    useradd -g genuser genuser && \
    mkdir /home/genuser && \
    chmod -R 777 /home/genuser

# copy dp_tools into container
COPY . /app

# set ownership and permissions
RUN chown -R genuser:genuser /app

# swith to user
USER genuser

RUN pip install /app &&  \
    echo "export PATH=/home/genuser/.local/bin:$PATH" >> ~/.bashrc

WORKDIR /home/genuser