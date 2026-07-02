# NucleoATAC — nucleosome calling from ATAC-seq (Schep et al., Genome Res 2015).
#
# Upstream GreenleafLab/NucleoATAC is Python 2.7-only and unmaintained. This
# image builds the community Python-3.11 modernization (GreenleafLab PR #99,
# alecethell/NucleoATAC@py311-modernization) pinned to its head commit, so the
# tool runs on a current Python/NumPy/SciPy/pysam stack.
#
# NOTE: pinned by commit SHA for reproducibility. amd64 is forced so the image
# also builds on Apple-Silicon hosts (emulated there), matching the other
# custom modules in this pipeline.
FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0

LABEL org.opencontainers.image.title="atacportray-nucleoatac"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.description="NucleoATAC nucleosome calling (Python 3.11 fork, GreenleafLab PR #99)"
LABEL org.opencontainers.image.licence="MIT"
LABEL org.opencontainers.image.source="https://github.com/alecethell/NucleoATAC/tree/py311-modernization"

# Runtime + build stack from conda (avoids compiling htslib for pysam). The
# Cython extensions (pyatac.fragments, nucleoatac.multinomial_cov) are compiled
# at install time against these exact numpy/cython versions.
RUN conda install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        "numpy>=1.26,<2.1" \
        "scipy>=1.11" \
        "pysam>=0.22" \
        "matplotlib>=3.8" \
        "cython>=3.0" \
        "setuptools>=68" \
        "samtools>=1.18" \
        wheel pip \
 && conda clean -afy

# Install NucleoATAC from the pinned py3.11 fork commit. --no-build-isolation
# reuses the conda-installed numpy/cython to build the extensions instead of
# pulling an isolated (possibly incompatible) build environment.
RUN pip install --no-cache-dir --no-build-isolation \
        "git+https://github.com/alecethell/NucleoATAC.git@22edcd0d026845c94013e12e5c550a17842edb7f"

ENV PATH=/opt/conda/bin:$PATH
# Headless matplotlib (NucleoATAC writes .eps diagnostic plots).
ENV MPLBACKEND=Agg
