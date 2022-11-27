FROM python:3.8.12
RUN python3 -m pip install pysam===0.18.0
RUN python3 -m pip install HTSeq===1.99.2
RUN cd /root && wget https://github.com/dieterich-lab/DCC/archive/refs/tags/v0.5.0.tar.gz && tar -xvzf v0.5.0.tar.gz && cd DCC-0.5.0 && python3 setup.py install
RUN python3 -m pip install pysam===0.18.0
RUN python3 -m pip install HTSeq===1.99.2
RUN apt-get update && apt-get install -y samtools


