#syntax=docker/dockerfile:1

FROM nvidia/cuda:11.6.1-runtime-ubuntu20.04

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y make git gcc g++

RUN git clone https://github.com/htgt/CRISPR-Analyser.git && \
    cd CRISPR-Analyser && \
    make

RUN CRISPR-Analyser/bin/crispr_analyser gather

COPY --chmod=0755 crispacuda .
