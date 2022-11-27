FROM staphb/fastqc:0.11.9
RUN apt-get update
RUN apt-get install -y flexbar
RUN apt-get install -y git
RUN cd /usr/local/bin/ && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && mv -vi bedtools.static.binary bedtools && chmod -v 755 bedtools
RUN apt-get install -y bzip2
RUN apt-get install -y gcc make g++ autoconf
RUN apt-get install -y libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev
RUN cd /root/ && wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && bunzip2  samtools-1.14.tar.bz2 && tar -xvf samtools-1.14.tar && cd samtools-1.14 && ./configure && make && make install
RUN apt-get install xz-utils  
RUN apt-get install -y libffi-dev
RUN apt-get install -y libreadline-dev
RUN apt-get install -y sqlite
RUN apt-get install -y libsqlite3-dev
RUN cd /root && wget  --no-check-certificate  https://www.openssl.org/source/openssl-1.1.1m.tar.gz && tar -xvzf openssl-1.1.1m.tar.gz  && cd openssl-1.1.1m && ./config && make && make install && ldconfig
RUN cd /root && wget https://www.python.org/ftp/python/3.10.1/Python-3.10.1.tar.xz && tar -xvf Python-3.10.1.tar.xz  &&  cd Python-3.10.1 && ./configure && make && make install



#bedtools (>= 2.25.0) and samtools (>= 1.3.1) 
#    1  wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
#    2  ls -alht
#    3  bunzip2
#    4  apt-get install -y bzip2
#    5  history 
