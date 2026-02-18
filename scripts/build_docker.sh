#!/bin/bash
#
# ElPGx Docker Build Script
# Builds and pushes all Docker images for production deployment
#
# Usage:
#   ./scripts/build_docker.sh           # Build only
#   ./scripts/build_docker.sh --push    # Build and push to Docker Hub
#

set -e

VERSION="1.1.0"
DOCKER_USER="elifffarslan"

echo "=============================================="
echo "ElPGx Docker Build Script v${VERSION}"
echo "=============================================="

cd "$(dirname "$0")/.."

# Check Docker is available
if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker is not installed."
    echo "Please install Docker Desktop: https://www.docker.com/products/docker-desktop"
    exit 1
fi

# Check Docker is running
if ! docker info &> /dev/null; then
    echo "ERROR: Docker daemon is not running."
    echo "Please start Docker Desktop."
    exit 1
fi

echo ""
echo "Building Docker images..."
echo ""

# Build ElPGx Python image
echo "[1/3] Building elpgx-python:${VERSION}..."
docker build \
    --platform linux/amd64 \
    -t ${DOCKER_USER}/elpgx-python:${VERSION} \
    -t ${DOCKER_USER}/elpgx-python:latest \
    -f docker/elpgx-python/Dockerfile \
    .

# Build ElPGx GATK image (optional - can use broadinstitute/gatk directly)
echo "[2/3] Building elpgx-gatk:${VERSION}..."
docker build \
    --platform linux/amd64 \
    -t ${DOCKER_USER}/elpgx-gatk:${VERSION} \
    -t ${DOCKER_USER}/elpgx-gatk:latest \
    -f docker/gatk/Dockerfile \
    .

# Build ElPGx HLA image
echo "[3/3] Building elpgx-hla:${VERSION}..."
docker build \
    --platform linux/amd64 \
    -t ${DOCKER_USER}/elpgx-hla:${VERSION} \
    -t ${DOCKER_USER}/elpgx-hla:latest \
    -f docker/hla/Dockerfile \
    .

echo ""
echo "=============================================="
echo "Docker images built successfully!"
echo "=============================================="
echo ""
docker images | grep elpgx
echo ""

# Push to Docker Hub if requested
if [[ "$1" == "--push" ]]; then
    echo "Pushing to Docker Hub..."
    echo ""

    # Login check
    if ! docker info 2>/dev/null | grep -q "Username"; then
        echo "Please login to Docker Hub first:"
        docker login
    fi

    echo "Pushing elpgx-python..."
    docker push ${DOCKER_USER}/elpgx-python:${VERSION}
    docker push ${DOCKER_USER}/elpgx-python:latest

    echo "Pushing elpgx-gatk..."
    docker push ${DOCKER_USER}/elpgx-gatk:${VERSION}
    docker push ${DOCKER_USER}/elpgx-gatk:latest

    echo "Pushing elpgx-hla..."
    docker push ${DOCKER_USER}/elpgx-hla:${VERSION}
    docker push ${DOCKER_USER}/elpgx-hla:latest

    echo ""
    echo "=============================================="
    echo "Docker images pushed to Docker Hub!"
    echo "=============================================="
    echo ""
    echo "Images available at:"
    echo "  - docker.io/${DOCKER_USER}/elpgx-python:${VERSION}"
    echo "  - docker.io/${DOCKER_USER}/elpgx-gatk:${VERSION}"
    echo "  - docker.io/${DOCKER_USER}/elpgx-hla:${VERSION}"
fi

echo ""
echo "To run the Nextflow pipeline:"
echo "  nextflow run main.nf \\"
echo "    --input <sample.bam> \\"
echo "    --reference <hg38.fa> \\"
echo "    --gene_tiers group1,group2,group3 \\"
echo "    -profile docker"
echo ""
