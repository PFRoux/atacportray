# atacportray/atacportray


[![GitHub Actions CI Status](https://github.com/atacportray/atacportray/actions/workflows/nf-test.yml/badge.svg)](https://github.com/atacportray/atacportray/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/atacportray/atacportray/actions/workflows/linting.yml/badge.svg)](https://github.com/atacportray/atacportray/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/atacportray/atacportray)

## Introduction

**atacportray/atacportray** is a bioinformatics pipeline for the multi-layer analysis of bulk **ATAC-seq** data from paired-end FASTQ files. Starting from raw reads and a reference genome, it produces a full epigenomic portrait (chromatin accessibility peaks, super-enhancers, coverage tracks and TF footprints) and, optionally, layers genetic and structural readouts mined from the same libraries: short-variant calls, copy-number profiles, telomere content and mitochondrial genotypes. All branches converge on a single [`MultiQC`](http://multiqc.info/) report.

The workflow is organised as a shared preprocessing core followed by toggleable analysis branches:

1. Read QC and adapter/quality trimming ([`fastp`](https://github.com/OpenGene/fastp))
2. **Epigenome branch** (`--run_epigenome`): Bowtie2 alignment, blacklist filtering, peak calling ([`MACS3`](https://github.com/macs3-project/MACS)), consensus peaks, coverage bigWigs ([`deepTools`](https://deeptools.readthedocs.io/)), super-enhancers ([`ROSE`](https://github.com/stjude/ROSE)) and optional TF footprinting ([`TOBIAS`](https://github.com/loosolab/TOBIAS))
3. **Variant branch** (`--run_variants`): BWA-MEM alignment, GATK MarkDuplicates + BQSR, multi-caller short-variant calling ([`DeepVariant`](https://github.com/google/deepvariant), [`FreeBayes`](https://github.com/freebayes/freebayes), [`GATK HaplotypeCaller`](https://gatk.broadinstitute.org/), [`bcftools`](https://samtools.github.io/bcftools/)), [`VEP`](https://www.ensembl.org/vep) annotation, [`vcf2maf`](https://github.com/mskcc/vcf2maf) conversion and a [`maftools`](https://github.com/PoisonAlien/maftools) oncoplot. This branch also yields the shared analysis-ready BAM consumed by the CNV / telomere / mito branches.
4. **CNV branch** (`--run_cnv`): copy-number profiling with [`QDNAseq`](https://bioconductor.org/packages/QDNAseq/)
5. **Telomere branch** (`--run_telomere`): telomere-content estimation with [`TelomereHunter`](https://www.dkfz.de/en/applied-bioinformatics/telomerehunter/telomerehunter.html)
6. **Mitochondrial branch** (`--run_mito`): mtDNA variant/heteroplasmy calling with [`mgatk`](https://github.com/caleblareau/mgatk)
7. **Regulatory network branch** (`--run_ananse`, `--run_coltron`): enhancer-based TF binding and GRN inference with [`ANANSE`](https://github.com/vanheeringen-lab/ANANSE), plus optional legacy [`Coltron`](https://pypi.org/project/coltron/) analysis from ROSE enhancer tables
8. Aggregate QC report ([`MultiQC`](http://multiqc.info/))

Only the epigenome branch runs by default; the genetic/structural and regulatory-network branches are opt-in because they require additional reference data and/or expression tables (and, for CNV/telomere/mito, `--run_variants` to build the shared analysis-ready BAM). See [usage documentation](docs/usage.md) for the branch dependency graph.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,condition,patient
IC_High_1,s1_R1.fastq.gz,s1_R2.fastq.gz,IC_High,P1
IC_Low_1,s2_R1.fastq.gz,s2_R2.fastq.gz,IC_Low,P2
```

Each row is one paired-end ATAC-seq sample. `condition` groups samples for consensus-peak building and downstream comparisons; `patient` is an optional grouping key. See [usage documentation](docs/usage.md) for the full column specification.

Now, you can run the pipeline using:

```bash
nextflow run atacportray/atacportray \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --fasta genome.fasta \
   --outdir <OUTDIR>
```

To enable the genetic/structural branches on a full-size reference:

```bash
nextflow run atacportray/atacportray \
   -profile <docker/singularity> \
   --input samplesheet.csv \
   --fasta genome.fasta \
   --run_variants --run_cnv --run_telomere --run_mito \
   --variant_callers deepvariant,freebayes,haplotypecaller,bcftools \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

atacportray/atacportray was originally written by Toumi et al..

We thank the following people for their extensive assistance in the development of this pipeline.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
