#!/bin/bash
nvcc -std=c++11 crispacuda.cu -o crispacuda && nvprof ./crispacuda -i /nfs/team87/crispr_indexes/GRCh38_index.bin -s 1106710986 -n 1000
