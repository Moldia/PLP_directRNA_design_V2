# Use Ubuntu 20.04 as base image
FROM ubuntu:20.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C

# Set working directory
WORKDIR /app

# Install system and bioinformatics build dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    make \
    gcc \
    g++ \
    python3 \
    python3-pip \
    build-essential \
    software-properties-common \
    automake \
    autoconf \
    perl \
    m4 \
    file \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    libdeflate-dev \
    bedtools \
    && apt-get clean

# Copy and install Python dependencies
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy and install local Python package
COPY PLP_directRNA_design_package /app/PLP_directRNA_design_package
WORKDIR /app/PLP_directRNA_design_package
RUN pip3 install .
WORKDIR /app  

# Download, build, and install samtools from source
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    && tar -xjf samtools-1.21.tar.bz2 \
    && rm samtools-1.21.tar.bz2 \
    && cd samtools-1.21 \
    && autoheader \
    && autoconf -Wno-syntax \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.21

# Default shell for debugging
CMD ["bash"]
