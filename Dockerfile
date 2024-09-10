#syntax=docker/dockerfile:1

FROM nvidia/cuda:11.6.1-runtime-ubuntu20.04

COPY --chmod=0755 crispacuda .
