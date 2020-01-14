# VSEARCH
#
# VERSION               2.4.3_sl

FROM      ubuntu
RUN apt-get update && apt-get install -y wget python3-pip
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip3 install pip --upgrade && pip3 install boto3 awscli
RUN mkdir /vsearch && mkdir /vsearch/src
RUN mkdir /fh
RUN mkdir /app
WORKDIR /vsearch/src
RUN wget https://github.com/torognes/vsearch/releases/download/v2.4.3/vsearch-2.4.3-linux-x86_64.tar.gz
RUN tar xzvf vsearch-2.4.3-linux-x86_64.tar.gz
RUN cp vsearch-2.4.3-linux-x86_64/bin/* /usr/local/bin/
RUN wget https://raw.githubusercontent.com/jgolob/sciluigi/containertask/tools/bucket_command_wrapper.py
RUN chmod +x bucket_command_wrapper.py && cp bucket_command_wrapper.py /usr/bin/
WORKDIR /vsearch
RUN rm -R /vsearch/src