# HLA Typing Image for PGx Pipeline
# Placeholder for future HLA module implementation (Group 4)

FROM python:3.11-slim

LABEL maintainer="PGx Genomize Team"
LABEL version="1.0.0"
LABEL description="HLA typing for PGx Pipeline (placeholder)"

# Install system dependencies for HLA typing
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libc-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Note: In production, this would install arcasHLA, OptiType, or HLA-HD
# Example for arcasHLA:
# RUN pip install arcasHLA

RUN pip install --no-cache-dir pysam>=0.22.0

# Placeholder message
RUN echo "HLA typing module is a placeholder. Implement arcasHLA/OptiType for production."

CMD ["echo", "HLA module placeholder - not yet implemented"]
