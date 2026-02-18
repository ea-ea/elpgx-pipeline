# ElPGx Pipeline Python Image
# Base image for Python-based pipeline stages
# Production-ready for Seven Bridges and other cloud platforms

FROM python:3.11-slim

LABEL maintainer="Elif Arslan"
LABEL version="1.1.0"
LABEL description="ElPGx - Clinical Pharmacogenetics Analysis Pipeline"
LABEL org.opencontainers.image.source="https://github.com/ea-ea/elpgx-pipeline"

# Install system dependencies + samtools + bcftools
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libc-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libffi-dev \
    make \
    samtools \
    bcftools \
    wget \
    unzip \
    openjdk-21-jre-headless \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install GATK
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip \
    && unzip gatk-4.5.0.0.zip -d /opt/ \
    && rm gatk-4.5.0.0.zip \
    && ln -s /opt/gatk-4.5.0.0/gatk /usr/local/bin/gatk

# Set working directory
WORKDIR /app

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    pysam>=0.22.0 \
    pydantic>=2.0.0 \
    cyvcf2>=0.30.0 \
    jsonschema>=4.0.0 \
    click>=8.0.0 \
    rich>=13.0.0

# Copy ElPGx package
COPY elpgx/ /app/elpgx/
COPY data/ /app/data/
COPY schemas/ /app/schemas/

# Set Python path
ENV PYTHONPATH=/app

# Verify installation
RUN python3 -c "from elpgx import __version__; print(f'ElPGx version: {__version__}')"

# No ENTRYPOINT — Nextflow needs a plain shell to run its scripts
CMD ["/bin/bash"]
