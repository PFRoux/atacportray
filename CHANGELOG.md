# atacportray/atacportray: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of atacportray/atacportray, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Shared preprocessing core: fastp trimming and reference indexing (Bowtie2 / BWA / faidx / GATK dict).
- **Epigenome branch** (`--run_epigenome`, default on): Bowtie2 alignment, blacklist filtering, MACS3 peak calling, per-condition consensus peaks, deepTools coverage bigWigs, ROSE super-enhancers and optional TOBIAS TF footprinting.
- **Variant branch** (`--run_variants`): BWA-MEM, GATK MarkDuplicates + BQSR, multi-caller short-variant calling (DeepVariant, FreeBayes, GATK HaplotypeCaller, bcftools), VEP annotation, vcf2maf conversion and a maftools oncoplot. Produces the shared analysis-ready BAM.
- **CNV branch** (`--run_cnv`): QDNAseq copy-number profiling.
- **Telomere branch** (`--run_telomere`): TelomereHunter content estimation.
- **Mitochondrial branch** (`--run_mito`): mgatk mtDNA variant calling.
- Custom local modules with containers, `meta.yml` and nf-test stubs: ROSE, TelomereHunter, mgatk, QDNAseq, TOBIAS (atacorrect/scorebigwig/bindetect), consensus_peaks, maftools_oncoplot.
- `test` (epigenome) and `test_full` (all branches) stub profiles.

### `Fixed`

### `Dependencies`

- Nextflow `>=24.10`
- nf-core/tools template 3.5.0

### `Deprecated`
