FROM jupyter/datascience-notebook:python-3.9.10
USER root
RUN apt-get update && apt-get install -y zip
COPY ct_aug.py /usr/local/bin/
CMD /bin/bash
