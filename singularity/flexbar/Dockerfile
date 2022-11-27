FROM ubuntu:20.04
RUN apt-get update
RUN apt-get install -y wget
RUN cd /root && wget https://github.com/seqan/flexbar/releases/download/v3.5.0/flexbar-3.5.0-linux.tar.gz && tar -xvzf flexbar-3.5.0-linux.tar.gz    &&  cd flexbar-3.5.0-linux && cp -vi ./flexbar /usr/local/bin && cp -vi *.so*   /usr/local/lib && ldconfig




