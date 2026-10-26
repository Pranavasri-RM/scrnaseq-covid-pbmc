FROM python:3.11-slim

# procps lets Nextflow measure CPU and memory use (same reason as the kb-python image)
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Reference data baked in, so tasks don't need internet at run time
RUN mkdir -p /opt/models /opt/genesets \
    && python -c "from celltypist import models; models.download_models(model=['Immune_All_Low.pkl'])" \
    && cp /root/.celltypist/data/models/Immune_All_Low.pkl /opt/models/ \
    && python -c "import gseapy as gp; lib = gp.get_library('MSigDB_Hallmark_2020', organism='Human'); open('/opt/genesets/MSigDB_Hallmark_2020.gmt', 'w').write(''.join(f'{t}\t\t' + '\t'.join(g) + '\n' for t, g in lib.items()))"

# Writable cache folders: Apptainer runs containers read-only, and numba and matplotlib cache on first use
ENV NUMBA_CACHE_DIR=/tmp/numba_cache \
    MPLCONFIGDIR=/tmp/matplotlib

LABEL org.opencontainers.image.source="https://github.com/Pranavasri-RM/scrnaseq-covid-pbmc"