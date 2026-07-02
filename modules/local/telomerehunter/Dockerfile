# TelomereHunter 1.1.0 (Feuerbach et al.) — python2 only, not on bioconda.
FROM continuumio/miniconda3:24.7.1-0

LABEL org.opencontainers.image.title="atacportray-telomerehunter"
LABEL org.opencontainers.image.version="1.1.0"
LABEL org.opencontainers.image.description="TelomereHunter telomere-content estimation (python2.7 + R plotting)"
LABEL org.opencontainers.image.licence="GPL-3.0"

RUN conda install -y -n base -c conda-forge -c bioconda \
        python=2.7 samtools=1.9 pysam r-base=3.6 r-ggplot2 r-reshape2 r-rcolorbrewer r-cowplot \
 && conda clean -afy

# TelomereHunter from PyPI (python2 wheel)
RUN pip install --no-cache-dir telomerehunter==1.1.0

ENV PATH=/opt/conda/bin:$PATH
