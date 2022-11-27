FROM rocker/rstudio@sha256:1271f20d89ff67fe55ed6be91632a76df96a405a9bcae1e8383781e464401ee4
RUN apt-get update
RUN apt-get install -y libxml2-dev  libfontconfig1-dev  libxt6  libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
RUN Rscript -e 'install.packages("devtools", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
RUN Rscript -e 'require(devtools);install_github("dieterich-lab/CircTest",ref="0f5a86c030dc35f5659fb3a7e9be2f16349790ff");library(CircTest)'
COPY circtest.cli.R /usr/local/bin/
COPY name_edit.py /usr/local/bin/
RUN chmod -v 755 /usr/local/bin/circtest.cli.R /usr/local/bin/name_edit.py
COPY circtest.cli.PLAIN.R /usr/local/bin/
RUN chmod -v 755 /usr/local/bin/circtest.cli.PLAIN.R
