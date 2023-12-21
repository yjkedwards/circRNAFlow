#!/bin/bash
n="kipoideeptarget"
en=${n}uab
docker build . -t local/${en}
sudo singularity build  --writable-tmpfs   ${en}.sif docker-daemon://local/${en}:latest

