#syntax=docker/dockerfile:1

FROM nvidia/cuda:11.6.1-runtime-ubuntu20.04

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y make git gcc g++ curl

RUN git clone https://github.com/htgt/CRISPR-Analyser.git && \
    cd CRISPR-Analyser && \
    make

RUN curl -vvv -O https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz & \
    gzip -d Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz & \
    CRISPR-Analyser/bin/crispr_analyser gather -i Homo_sapiens.GRCh38.dna.chromosome.1.fa -o chrom_1.tsv -e 1

COPY --chmod=0755 crispacuda .
