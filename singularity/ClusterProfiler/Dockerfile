FROM r-base:4.2.1
#FROM jupyter/datascience-notebook:r-4.2
USER root
RUN apt-get update
RUN apt-get install -y libfreetype-dev
RUN apt-get install -y libfreetype6 
RUN apt-get install -y libz-dev
RUN apt-get install -y libxml2-dev
RUN apt-get install -y libglpk-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y curl
COPY install_cp.R /root/
#RUN cd /root && Rscript install_cp.R
COPY install_o.R /root/
RUN cd /root && Rscript install_o.R
RUN apt-get install -y libxtst6 libxt6
RUN apt-get install -y libcurl4-openssl-dev
RUN cd /root && Rscript install_o.R
#install GNU parallel
RUN cd /root/ && wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2 && tar -xvf parallel-latest.tar.bz2 && cd `ls -1t|grep -Pi 'para'|tail --lines=1` && ./configure && make && make install
#freetype2 header files ended up not being found, this cp command seems to fix it
#RUN apt-get install -y libfreetype-dev libfreetype6 libfreetype6-dev
#RUN cd /usr/include/freetype2/ && cp -vr * /opt/conda/include/
RUN apt-get install -y libfontconfig1-dev
RUN cd /root/ && Rscript install_cp.R
RUN apt-get install -y python3-dev
RUN apt-get install -y pip
RUN python3 -m pip   install pandas  --break-system-packages
COPY install_cp.4.4.R /root/
RUN cd /root && Rscript install_cp.4.4.R
COPY cp_script.R /usr/local/bin/
CMD /bin/bash
