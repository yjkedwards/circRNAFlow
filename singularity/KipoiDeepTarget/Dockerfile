FROM kipoi/kipoi-docker:sharedpy3keras2tf1
RUN useradd -ms /bin/bash kipoi_user
USER kipoi_user
ENV HOME=/home/kipoi_user
RUN python3 -m pip install biopython

#load in deeptarget model/files
RUN kipoi env create deepTarget
RUN source activate kipoi-deepTarget && python3 -c "import kipoi;model = kipoi.get_model('deepTarget');"

#attempt to disable calls/netIO with github (don't overwrite models/python files etc)
RUN find /home  -type f -iname "*.py" | xargs -tI {} sed -r -i "s/github\.com/localhost/g" {}
RUN find /home  -type f -iname "*.yaml" | xargs -tI {} sed -r -i "s/github\.com/localhost/g" {}
USER root
RUN find /opt -type f -iname "*.py" | xargs -tI {} sed -r -i "s/github\.com/localhost/g" {}
RUN find /opt -type f -iname "*.yaml" | xargs -tI {} sed -r -i "s/github\.com/localhost/g" {}
USER kipoi_user
RUN sed -ir -r "s/github\.com/localhost/g"   /home/kipoi_user/.kipoi/models/.git/config
COPY deep_target_runner.py /usr/local/bin/

#export HOME as /home/kipoi_user
