#syntax=docker/dockerfile:1

FROM nvidia/cuda:11.6.1-runtime-ubuntu20.04

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y make git gcc g++ curl

RUN git clone https://github.com/htgt/CRISPR-Analyser.git && \
    cd CRISPR-Analyser && \
    make

RUN curl -vvv -O https://ftp.sanger.ac.uk/pub/teams/229/crispr_indexes/GRCh38_index.bin.gz && \
    gzip -d GRCh38_index.bin.gz

COPY --chmod=0755 crispacuda .

CMD ["./crispacuda", "-i", "GRCh38_index.bin", "-p", "8000"]
