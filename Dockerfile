#
# PPaxe standalone docker requires
#    + debian:jessie or ubuntu:16.04
#    + python 2.7
#
# Build this docker with:
#   docker build -t ppaxe.docker -f=./Dockerfile_ppaxe .
#
# Run this docker with:
#   docker run -v /local/path/to/output:/ppaxe/output:rw \
#               ppaxe.docker -v -p ./papers.pmids -o ./papers.tbl -r ./papers.html
#   the container working local folder is set to /ppaxe/output
#   where the program will return by default the results.
#   You must mount the container folder to your "/local/path/to/output"
#   to keep the final results, with the docker "-v" switch.
#
FROM ubuntu:16.04

MAINTAINER Josep F Abril, jabril@ub.edu
MAINTAINER Sergio Castillo-Lara, s.cast.lara@gmail.com

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
               build-essential \
               ca-certificates \
               gcc \
               git \
               libpq-dev \
               make \
               python-pip \
               python2.7 \
               python2.7-dev \
               wget \
               ssh \
               unzip \
               openjdk-8-jdk \
    && apt-get autoremove \
    && apt-get clean

RUN pip install -U "pycorenlp==0.3.0"
RUN pip install -U "scipy==0.17.0"
RUN pip install -U "sklearn==0.0"
RUN pip install -U "requests==2.4.3"
RUN pip install -U "scikit-learn==0.18.2"
RUN pip install -U "matplotlib==2.0.2"

RUN wget http://nlp.stanford.edu/software/stanford-corenlp-full-2017-06-09.zip \
    && unzip stanford-corenlp-full-2017-06-09.zip \
    && wget https://www.dropbox.com/s/ec3a4ey7s0k6qgy/FINAL-ner-model.AImed%2BMedTag%2BBioInfer.ser.gz?dl=0 -O FINAL-ner-model.AImed+MedTag+BioInfer.ser.gz \
    && wget http://nlp.stanford.edu/software/stanford-english-corenlp-2017-06-09-models.jar \
         -O stanford-corenlp-full-2017-06-09/stanford-english-corenlp-2017-06-09-models.jar 

RUN git clone https://github.com/CompGenLabUB/ppaxe.git

WORKDIR /ppaxe

RUN wget https://www.dropbox.com/s/t6qcl19g536c0zu/RF_scikit.pkl?dl=0 \
      -O ./ppaxe/data/RF_scikit.pkl \
    && pip install ./

# WORKDIR /stanford-corenlp-full-2017-06-09
# 
# RUN java -mx1000m \
#          -cp ./stanford-corenlp-3.8.0.jar:stanford-english-corenlp-2017-06-09-models.jar \
#          edu.stanford.nlp.pipeline.StanfordCoreNLPServer \
#          -port 9000 \
#          -serverProperties ../ppaxe/ppaxe/data/server.properties &

RUN mkdir -vp /ppaxe/output \
    && chmod -v a+rwx /ppaxe/output

RUN \
  echo '#!/bin/bash\n\
SPD="/stanford-corenlp-full-2017-06-09"\n\
java -mx1000m \\n\
         -cp $SPD/stanford-corenlp-3.8.0.jar:stanford-english-corenlp-2017-06-09-models.jar \\n\
         $SPD/edu.stanford.nlp.pipeline.StanfordCoreNLPServer \\n\
         -port 9000 \\n\
         -serverProperties /ppaxe/ppaxe/data/server.properties &\n\
\n\
cd /ppaxe/output\n\
/ppaxe/bin/ppaxe $@\n' \
  > /ppaxe/entrypoint.sh \
  && chmod +x /ppaxe/entrypoint.sh

WORKDIR /ppaxe/output

ENTRYPOINT ["/ppaxe/entrypoint.sh"]
