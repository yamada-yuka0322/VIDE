FROM debian:stable


RUN apt update && apt install -y git gcc g++ python3 python3-dev python3-pip python3-numpy python3-setuptools cmake libz-dev coreutils m4

WORKDIR /build

###ADD vide_key_ssh .
ENV GIT_SSH_COMMAND='ssh -i /ssh-key/vide_key_ssh -o "StrictHostKeyChecking no"'
RUN pip3 install --no-cache-dir 'astropy<4.1'
RUN git clone --recursive git@bitbucket.org:glavaux/vide_next.git 

WORKDIR /build/vide_next

RUN git checkout python3
#RUN python3 setup.py build
RUN pip3 install --verbose .
