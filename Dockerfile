#
# PPaxe standalone docker requires
#    + debian:jessie or ubuntu:16.04
#    + python 2.7
#
# Build this docker with:
#   docker build -t ppaxe.docker -f=./Dockerfile .
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
               apt-utils \
               python-pip \
               python2.7 \
               python2.7-dev \
               wget \
               ssh \
               unzip \
               openjdk-8-jdk \
    && apt-get autoremove \
    && apt-get clean

RUN pip install --upgrade pip
RUN pip install -U "pycorenlp==0.3.0"
RUN pip install -U "scipy==0.17.0"
RUN pip install -U "sklearn==0.0"
RUN pip install -U "requests==2.4.3"
RUN pip install -U "scikit-learn==0.18.2"
RUN pip install -U "matplotlib==2.0.2"
RUN pip install -U "networkx==2.2"
# RUN pip install -U "pysqlite"  # this is only for webapp
RUN pip install -U "uuid"

RUN wget http://nlp.stanford.edu/software/stanford-corenlp-full-2017-06-09.zip \
    && unzip stanford-corenlp-full-2017-06-09.zip \
    && wget https://compgen.bio.ub.edu/datasets/PPaxe_files/FINAL-ner-model.AImed%2BMedTag%2BBioInfer.ser.gz \
         -O /stanford-corenlp-full-2017-06-09/FINAL-ner-model.AImed+MedTag+BioInfer.ser.gz \
    && wget http://nlp.stanford.edu/software/stanford-english-corenlp-2017-06-09-models.jar \
         -O /stanford-corenlp-full-2017-06-09/stanford-english-corenlp-2017-06-09-models.jar 

RUN echo "# Installing ppaxe lib..." \
    && git --no-pager clone https://github.com/scastlara/ppaxe.git
#        # repo clone at https://github.com/CompGenLabUB/ppaxe.git
RUN sed -i 's%\.\./%/stanford-corenlp-full-2017-06-09/%' \
           /ppaxe/ppaxe/data/server.properties

WORKDIR /ppaxe

RUN wget https://compgen.bio.ub.edu/datasets/PPaxe_files/RF_scikit.pkl \
      -O ./ppaxe/data/RF_scikit.pkl \
    && pip install ./

RUN mkdir -vp /ppaxe/output \
    && chmod -v a+rwx /ppaxe/output

##
## System params... set up your own \n\
##
ENV CORENLP_THREADS=4
ENV CORENLP_MAXMEM=4g

RUN \
  echo '#!/bin/bash\n\
SPD="/stanford-corenlp-full-2017-06-09"\n\
nohup java -Xms${CORENLP_MAXMEM} -Xmx${CORENLP_MAXMEM} \\\n\
         -cp $SPD/stanford-corenlp-3.8.0.jar:$SPD/stanford-english-corenlp-2017-06-09-models.jar \\\n\
         edu.stanford.nlp.pipeline.StanfordCoreNLPServer \\\n\
         -port 9000 -threads ${CORENLP_THREADS} \\\n\
         -serverProperties /ppaxe/ppaxe/data/server.properties \\\n\
         2> /dev/null 1>&2 &\n\
\n\
cd /ppaxe/output\n\
/ppaxe/bin/ppaxe $@\n' \
  > /ppaxe/entrypoint.sh \
  && chmod +x /ppaxe/entrypoint.sh

# RUN cat /ppaxe/entrypoint.sh

WORKDIR /ppaxe/output

ENTRYPOINT ["/ppaxe/entrypoint.sh"]
