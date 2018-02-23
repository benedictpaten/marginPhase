FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# update and install dependencies
RUN apt-get update && apt-get -y install git make gcc bzip2 lzma-dev zlib1g-dev libcurl4-openssl-dev libcrypto++-dev \
    libpthread-stubs0-dev g++ wget libbz2-dev liblzma-dev autoconf unzip php5-mcrypt python-pip python-dev
RUN pip install numpy pysam

# build and install cPecan
WORKDIR /opt
ENV ROOT /opt
ENV CPECAN_COMMIT 5acfa1eb37569c828b5ecd0f8f812510b713c770
ENV SONLIB_COMMIT 265d6c510cff90412d88e5ab9945818cfdff6e41
RUN wget https://github.com/benedictpaten/sonLib/archive/$SONLIB_COMMIT.zip && unzip $SONLIB_COMMIT.zip \
    && cd $ROOT/sonLib-$SONLIB_COMMIT && make && cd $ROOT && mv $ROOT/sonLib-$SONLIB_COMMIT $ROOT/sonLib \
    && wget https://github.com/benedictpaten/cPecan/archive/$CPECAN_COMMIT.zip && unzip $CPECAN_COMMIT.zip \
    && cd $ROOT/cPecan-$CPECAN_COMMIT && make && cd $ROOT \
    && mv $ROOT/sonLib $ROOT/sonLib-$SONLIB_COMMIT \
    && mkdir -p $ROOT/cPecan && mv $ROOT/sonLib-$SONLIB_COMMIT/bin/cPecan* $ROOT/cPecan \
    && mv $ROOT/cPecan-$CPECAN_COMMIT/*.py $ROOT/cPecan \
    && rm -rf $ROOT/*$SONLIB_COMMIT* $ROOT/*$CPECAN_COMMIT*

# install latest cmake
WORKDIR /tmp
RUN mkdir /opt/cmake && \
    wget https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh && \
    sh /tmp/cmake-3.7.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    mv /opt/cmake/bin/cmake /usr/local/bin/cmake && \
    rm -rf /opt/cmake

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
