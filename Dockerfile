# Use Ubuntu as the base image
FROM ubuntu:20.04

# Set environment variables to prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    make \
    gcc \
    python3 \
    python3-pip \
    python3-cutadapt \
    && apt-get clean

# Install Python dependencies
WORKDIR /app
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Install bedtools
RUN apt-get update && apt-get install -y \
    software-properties-common \
    && add-apt-repository universe \
    && apt-get update \
    && apt-get install -y bedtools \
    && apt-get clean

# Add local Python package (PLP_directRNA_design)
COPY PLP_directRNA_design /app/PLP_directRNA_design

# Expose the port for Jupyter Notebook
EXPOSE 8888

# Install Jupyter Notebook
RUN pip3 install notebook

# Set the working directory
WORKDIR /app

# Start a shell for manual debugging
CMD ["bash"]

