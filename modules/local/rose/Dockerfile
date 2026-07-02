# ROSE — Rank Ordering of Super-Enhancers (Young Lab)
# Not available on bioconda; built here from source for atacportray.
FROM continuumio/miniconda3:24.7.1-0

LABEL org.opencontainers.image.title="atacportray-rose"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.description="ROSE super-enhancer caller (python2.7 + samtools + bamliquidator)"
LABEL org.opencontainers.image.licence="MIT"

RUN conda config --add channels defaults \
 && conda config --add channels bioconda \
 && conda config --add channels conda-forge \
 && conda install -y -n base -c conda-forge -c bioconda \
        python=2.7 samtools=1.9 r-base=3.6 bamliquidator numpy \
 && conda clean -afy

# Fetch ROSE from the maintained mirror and expose the scripts on PATH
RUN mkdir -p /opt/rose \
 && cd /opt/rose \
 && wget -q https://bitbucket.org/young_computation/rose/get/HEAD.tar.gz -O rose.tar.gz \
 && tar --strip-components=1 -xzf rose.tar.gz \
 && rm rose.tar.gz \
 && chmod +x *.py

ENV PATH=/opt/rose:/opt/conda/bin:$PATH

# Helper: convert narrowPeak/BED to ROSE GFF (shipped via module resources)
COPY resources/usr/bin/rose_peaks_to_gff.py /usr/local/bin/rose_peaks_to_gff.py
RUN chmod +x /usr/local/bin/rose_peaks_to_gff.py
