#!/bin/bash
n="star"
en=${n}uab
docker build . -t local/${en}
sudo singularity build ${en}.sif docker-daemon://local/${en}:latest

