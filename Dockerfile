FROM nvcr.io/nvidia/cuda:11.7.1-devel-ubuntu20.04

COPY imagebuild/sources.list /etc/apt/sources.list

# Install prerequisites
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get -y install \
    autoconf \
    automake \
    build-essential \
    gfortran \
    git \
    python2.7 \
    python-dev \
    python3-dev \
    latex2html \
    libcfitsio-bin \
    libcfitsio-dev \
    libfftw3-bin \
    libfftw3-dev \
    libglib2.0-dev \
    libpng-dev \
    libtool \
    libx11-dev \
    pgplot5 \
    tcsh \
    libbz2-dev \
    python-tk \
    wget && \
    apt-get clean all && \
    rm -r /var/lib/apt/lists/*

# Add pgplot environment variables
ENV PGPLOT_DIR=/usr/local/pgplot
ENV PGPLOT_DEV=/Xserve

ADD . /home/soft/presto
ENV PRESTO /home/soft/presto
ENV LD_LIBRARY_PATH /home/soft/presto/lib

WORKDIR /home/soft/presto

RUN python imagebuild/get-pip2.py && python3 imagebuild/get-pip.py

RUN pip2 install numpy && \
    pip2 install minio \
    pyfits \
    fitsio \
    matplotlib \
    scipy \
    astropy \
    future

# Install python dependancies
RUN pip3 install numpy \
    minio \
    matplotlib \
    scipy \
    astropy \
    pymongo \
    boto3 \
    future 

# Install presto python scripts

WORKDIR /home/soft/presto/src
# The following is necessary if your system isn't Ubuntu 20.04
RUN make cleaner
# Now build from scratch
RUN make libpresto
WORKDIR /home/soft/presto
RUN pip3 install /home/soft/presto && \
    sed -i 's/env python/env python3/' /home/soft/presto/bin/*py && \
    python3 tests/test_presto_python.py


# Installs all the C dependancies -----------------------------
WORKDIR /home/soft

# Install psrcat
RUN wget https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz && \
    gunzip psrcat_pkg.tar.gz && \
    tar -xvf psrcat_pkg.tar && \
    rm psrcat_pkg.tar && \
    cd psrcat_tar && \
    ls && \
    bash makeit && \
    cp psrcat /usr/bin
ENV PSRCAT_FILE /home/soft/psrcat_tar/psrcat.db

# Install tempo
# RUN git clone https://github.com/nanograv/tempo.git && \
#     cd tempo && \
#     ./prepare && \
#     ./configure && \
#     make && \
#     make install

# Install tempo
COPY imagebuild/tempo /home/soft/tempo

RUN cd tempo && \
    ./prepare && \
    ./configure && \
    make && \
    make install

ENV TEMPO /home/soft/tempo


WORKDIR /home/soft/presto/src
# Build presto
RUN make makewisdom && \
    make clean && \
    make binaries

WORKDIR /home/soft/
ENV PATH="/home/soft/presto/bin/:${PATH}"

