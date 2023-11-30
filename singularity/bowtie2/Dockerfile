FROM python:3.8.12
#FROM ubuntu:20.04
RUN apt-get update
RUN apt-get install -y wget unzip 
RUN apt-get install -y libtbb-dev
COPY bowtie2-2.3.0-linux-x86_64.zip /root/
RUN cd /root && unzip bowtie2-2.3.0-linux-x86_64.zip
RUN cd /root && mv -v bowtie2-2.3.0 /usr/local && cd /usr/local/bin &&  find ../bowtie2-2.3.0/ -maxdepth 1   -type f -iname "bow*" | xargs -tI {} ln -vs {}
COPY bowtie2 /usr/local/bowtie2-2.3.0/bowtie2
CMD /bin/bash
