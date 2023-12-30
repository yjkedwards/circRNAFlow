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
RUN chmod -vR 777 /home/kipoi_user/.theano/
RUN find /opt -type f -iname "*.py" | xargs -tI {} sed -r -i "s/github\.com/localhost/g" {}
RUN find /opt -type f -iname "*.yaml" | xargs -tI {} sed -r -i "s/github\.com/localhost/g" {}


#attempt to set home as /home/kipoi_user instead of using ~
#root@e662729f0e7a:/opt/conda# find .|grep -Pi 'kipoi/config.py'
#./envs/kipoi-shared__envs__kipoi-py3-keras2-tf1/lib/python3.7/site-packages/kipoi/config.py
#./envs/kipoi-env/lib/python3.10/site-packages/kipoi/config.py
COPY config.py /opt/conda/envs/kipoi-shared__envs__kipoi-py3-keras2-tf1/lib/python3.7/site-packages/kipoi/config.py
COPY config.py /opt/conda/envs/kipoi-env/lib/python3.10/site-packages/kipoi/config.py
COPY config.py /home/kipoi_user/.conda/envs/kipoi-deepTarget/lib/python3.7/site-packages/kipoi/config.py
COPY sources.py /opt/conda/envs/kipoi-env/lib/python3.10/site-packages/kipoi/sources.py
COPY sources.py /opt/conda/envs/kipoi-shared__envs__kipoi-py3-keras2-tf1/lib/python3.7/site-packages/kipoi/sources.py
COPY sources.py /home/kipoi_user/.conda/envs/kipoi-deepTarget/lib/python3.7/site-packages/kipoi/sources.py
COPY model.py /opt/conda/envs/kipoi-shared__envs__kipoi-py3-keras2-tf1/lib/python3.7/site-packages/kipoi/model.py
COPY model.py /opt/conda/envs/kipoi-env/lib/python3.10/site-packages/kipoi/model.py
COPY model.py /home/kipoi_user/.conda/envs/kipoi-deepTarget/lib/python3.7/site-packages/kipoi/model.py
USER kipoi_user


RUN sed -ir -r "s/github\.com/localhost/g"   /home/kipoi_user/.kipoi/models/.git/config
COPY deep_target_runner.py /usr/local/bin/
COPY m_patch.txt /home/kipoi_user/
RUN cat /home/kipoi_user/m_patch.txt >> /home/kipoi_user/.kipoi/models/deepTarget/bio_utils.py
USER root
RUN find /home/kipoi_user/.* -maxdepth 0 -type d|grep -Piv '\.$'|xargs -tI {} chmod -vR 777 {}
RUN chmod -vR 777 /home/kipoi_user/
#make all of the user's file owned by root to avoid errors/warnings about unknown user ID
RUN chown -vR root:root /home/kipoi_user/
USER kipoi_user
