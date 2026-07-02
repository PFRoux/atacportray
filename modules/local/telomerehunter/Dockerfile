# TelomereHunter 1.1.0 (Feuerbach et al.) is python2 only and not on bioconda.
FROM python:2.7.18-buster

LABEL org.opencontainers.image.title="atacportray-telomerehunter"
LABEL org.opencontainers.image.version="1.1.0"
LABEL org.opencontainers.image.description="TelomereHunter telomere-content estimation"
LABEL org.opencontainers.image.licenses="GPL-3.0"
LABEL org.opencontainers.image.source="https://github.com/PFRoux/atacportray"

RUN sed -i 's/deb.debian.org/archive.debian.org/g; s/security.debian.org/archive.debian.org/g; /buster-updates/d' /etc/apt/sources.list \
 && apt-get -o Acquire::Check-Valid-Until=false update \
 && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        samtools \
        zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        numpy==1.16.6 \
        pysam==0.15.4 \
        PyPDF2==1.26.0 \
        telomerehunter==1.1.0
