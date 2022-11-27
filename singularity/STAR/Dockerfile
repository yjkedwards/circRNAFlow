FROM ubuntu:20.04
RUN apt-get update
RUN apt-get install -y wget
RUN cd /root && wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz && tar -xvzf 2.7.9a.tar.gz &&  cd /root/STAR-2.7.9a/bin/Linux_x86_64_static && mv -v STAR* /usr/local/bin/

