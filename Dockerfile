FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    bowtie2 \
    openjdk-8-jdk \
    pbzip2 \
    pigz \
    python3.10 \
    python3-pip \
    python3-biopython \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install RSeQC
RUN pip3 install --no-cache-dir RSeQC

# Install Atria
RUN wget -q https://github.com/cihga39871/Atria/releases/download/v3.1.2/atria-3.1.2-linux.tar.gz \
    && tar xf atria-3.1.2-linux.tar.gz \
    && mv atria-3.1.2 /opt/atria

# Install findtail
RUN wget -q https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/findtail/findtail_v1.01 \
    && chmod +x findtail_v1.01 \
    && mv findtail_v1.01 /usr/local/bin/findtail_v1.01

# Add the Atria executable to the PATH
ENV PATH="/opt/atria/bin:${PATH}"