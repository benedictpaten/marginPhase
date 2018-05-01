FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# update and install dependencies
RUN apt-get update && apt-get -y install git make wget autoconf gcc g++ bzip2 libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev liblzma-dev libssl-dev

# install latest cmake
WORKDIR /tmp
RUN wget https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh
RUN mkdir /opt/cmake
RUN sh /tmp/cmake-3.7.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license
RUN ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

# get marginPhase
WORKDIR /opt
ADD tempMarginPhase /opt/marginPhase

# build htslib
WORKDIR /opt/marginPhase/externalTools/htslib
RUN autoconf ; autoheader ; ./configure ; make clean ; make

# build sonlib
WORKDIR /opt/marginPhase/externalTools/sonLib
RUN make clean ; make

# build marginPhase
WORKDIR /opt/marginPhase
RUN cmake .
RUN make clean && make

# setup entrypoint
COPY wrapper.sh /opt/
ENTRYPOINT ["sh", "/opt/wrapper.sh"]
