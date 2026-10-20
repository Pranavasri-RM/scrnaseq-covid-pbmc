# Base: official slim Python image (Debian Linux with Python 3.12)
FROM python:3.12-slim

# System libraries kallisto/bustools need, plus procps for Nextflow
RUN apt-get update \
 && apt-get install -y --no-install-recommends libstdc++6 zlib1g procps \
 && rm -rf /var/lib/apt/lists/*

# kb-python, pinned to an exact version
RUN pip install --no-cache-dir kb-python==0.30.2

# Link this image to the GitHub repo (shows under the repo's Packages)
LABEL org.opencontainers.image.source="https://github.com/Pranavasri-RM/scrnaseq-covid-pbmc"