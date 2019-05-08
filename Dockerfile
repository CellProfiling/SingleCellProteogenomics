FROM conda/miniconda3
LABEL maintainer="Anthony Cesnik <cesnik@wisc.edu>"

WORKDIR /app
COPY . ./
RUN apt-get update && \
    pip install scanpy pandas && \
    conda init && \
    conda install -y python-louvain