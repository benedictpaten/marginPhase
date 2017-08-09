FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# update and install dependencies
RUN apt-get update && apt-get -y install git make gcc bzip2 lzma-dev zlib1g-dev libcurl4-openssl-dev libcrypto++-dev libpthread-stubs0-dev g++ wget libbz2-dev liblzma-dev autoconf

# install latest cmake
WORKDIR /tmp
RUN wget https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh
RUN mkdir /opt/cmake
RUN sh /tmp/cmake-3.7.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license
RUN ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

# get marginPhase
WORKDIR /opt
#TODO don't clone it, copy all files from current toil directory?
RUN git clone  --recursive https://github.com/benedictpaten/marginPhase.git

# build htslib
WORKDIR /opt/marginPhase/externalTools/htslib
RUN autoconf ; autoheader ; ./configure ; make

# build marginPhase
WORKDIR /opt/marginPhase
RUN cmake .
RUN make

# setup entrypoint
COPY wrapper.sh /opt/
ENTRYPOINT ["sh", "/opt/wrapper.sh"]
