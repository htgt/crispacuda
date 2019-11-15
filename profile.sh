#!/bin/bash
nvcc -std=c++11 crispacuda.cu -o crispacuda && nvprof ./crispacuda
