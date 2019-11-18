#!/bin/bash
nvcc -std=c++11 crispacuda.cu -o crispacuda && cat ~/268735.txt | xargs nvprof ./crispacuda -i /nfs/team87/crispr_indexes/GRCh38_index.bin
#./crispacuda -i /nfs/team87/crispr_indexes/GRCh38_index.bin -s 1106711098 -n 500
