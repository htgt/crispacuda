#!/bin/bash
nvcc -std=c++11 main.cu -o crispacuda && nvprof ./crispacuda
