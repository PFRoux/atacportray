# atacportray/atacportray: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [fastp](#fastp) - Adapter/quality trimming and read QC
- [Epigenome branch](#epigenome-branch)
  - [Bowtie2](#bowtie2) - Read alignment
  - [MACS3](#macs3) - Peak calling
  - [Consensus peaks](#consensus-peaks) - Per-condition reproducible peak sets
  - [deepTools](#deeptools) - Coverage bigWig tracks
  - [ROSE](#rose) - Super-enhancer identification
  - [NucleoATAC](#nucleoatac) - Nucleosome positioning (optional, `--run_nucleoatac`)
  - [TOBIAS](#tobias) - TF footprinting (optional, `--run_footprinting`)
- [Variant branch](#variant-branch) (`--run_variants`)
  - [BWA-MEM + GATK](#bwa-mem--gatk) - Alignment, MarkDuplicates, BQSR (analysis-ready BAM)
  - [Variant callers](#variant-callers) - DeepVariant, FreeBayes, HaplotypeCaller, bcftools
  - [VEP + vcf2maf + maftools](#annotation) - Annotation, MAF conversion, oncoplot
- [CNV branch](#cnv-branch) (`--run_cnv`) - QDNAseq copy-number profiles
- [Telomere branch](#telomere-branch) (`--run_telomere`) - TelomereHunter content estimates
- [Mitochondrial branch](#mitochondrial-branch) (`--run_mito`) - mgatk mtDNA variants
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Preprocessing

#### fastp

<details markdown="1">
<summary>Output files</summary>

- `preprocessing/fastp/`
  - `*.fastp.html`, `*.fastp.json` - Per-sample trimming reports
  - `*.fastp.log` - Run log

</details>

[fastp](https://github.com/OpenGene/fastp) trims adapters and low-quality bases and reports read-level QC, aggregated in MultiQC.

### Epigenome branch

#### Bowtie2

<details markdown="1">
<summary>Output files</summary>

- `epigenome/bowtie2/`
  - `*.bam`, `*.bai` - Filtered, sorted, blacklist-cleaned alignments

</details>

#### MACS3

<details markdown="1">
<summary>Output files</summary>

- `epigenome/macs3/`
  - `*_peaks.narrowPeak`, `*_peaks.xls`, `*_summits.bed` - Called accessibility peaks

</details>

#### Consensus peaks

<details markdown="1">
<summary>Output files</summary>

- `epigenome/consensus/`
  - `*.consensus.bed` - Reproducible peak set per `condition`

</details>

#### deepTools

<details markdown="1">
<summary>Output files</summary>

- `epigenome/coverage/`
  - `*.bw` - CPM-normalised coverage tracks

</details>

#### ROSE

<details markdown="1">
<summary>Output files</summary>

- `epigenome/rose/`
  - `*_SuperStitched.table.txt` - Ranked super-enhancers (stitched regions above the cutoff)
  - `*_AllStitched.table.txt` - All stitched enhancers with signal ranking and super-enhancer status
  - `*_Stitched_withSuper.bed` - BED of stitched enhancers flagged as super-enhancers
  - `*_SuperStitched.table_withGENES.txt` - Super-enhancers with nearby genes assigned by ROSE_geneMapper
  - `*_Plot_points.png` - Hockey-stick ranking plot

</details>

#### NucleoATAC

<details markdown="1">
<summary>Output files</summary>

- `epigenome/nucleoatac/<sample>/`
  - `*.nucpos.bed.gz` - Called nucleosome dyad positions
  - `*.nucmap_combined.bed.gz` - Combined nucleosome map
  - `*.nfrpos.bed.gz` - Nucleosome-free regions
  - `*.occpeaks.bed.gz` - Occupancy peak calls
  - `*.occ.bedgraph.gz` - Nucleosome occupancy track
  - `*.nucleoatac_signal.bedgraph.gz`, `*.nucleoatac_signal.smooth.bedgraph.gz` - NucleoATAC signal tracks

</details>

Produced only when `--run_nucleoatac` is set. NucleoATAC calls nucleosome
positions and nucleosome-free regions within the per-sample MACS3 peak regions
using the analysis-ready BAM and the reference FASTA (for Tn5-bias correction).
Runs in a Python 3.11 container built from the community fork (see
[containers.md](containers.md)).

#### TOBIAS

<details markdown="1">
<summary>Output files</summary>

- `epigenome/tobias/`
  - `*_corrected.bw`, `*_footprints.bw` - Tn5-bias-corrected signal and footprint scores
  - `bindetect/` - Per-motif differential binding tables and plots

</details>

Produced only when `--run_footprinting` is set and `--tobias_motifs` is provided.

### Variant branch

#### BWA-MEM + GATK

<details markdown="1">
<summary>Output files</summary>

- `variants/bam/<sample>/`
  - `*.recalibrated.bam`, `*.bai` - Analysis-ready BAM shared with CNV/telomere/mito branches

</details>

#### Variant callers

<details markdown="1">
<summary>Output files</summary>

- `variants/calls/<caller>/`
  - `*.vcf.gz`, `*.vcf.gz.tbi` - Per-caller short-variant calls (DeepVariant, FreeBayes, HaplotypeCaller, bcftools)

</details>

#### Annotation

<details markdown="1">
<summary>Output files</summary>

- `variants/annotated/`
  - `*.vep.vcf.gz` - VEP-annotated variants
  - `*.maf` - MAF conversion (vcf2maf)
  - `oncoplot.pdf` - maftools oncoplot across samples (when `--run_oncoplot`)

</details>

### CNV branch

<details markdown="1">
<summary>Output files</summary>

- `cnv/qdnaseq/`
  - `*_copynumber.igv`, `*_calls.txt` - Binned log2 ratios and gain/loss calls
  - `*_profile.png` - Genome-wide copy-number plot

</details>

### Telomere branch

<details markdown="1">
<summary>Output files</summary>

- `telomere/telomerehunter/<sample>/`
  - `*_summary.tsv` - Telomere-content estimates (intratelomeric reads per million)
  - `*_plots/` - Telomere-content and pattern plots

</details>

### Mitochondrial branch

<details markdown="1">
<summary>Output files</summary>

- `mito/mgatk/<sample>/`
  - `*.variant_stats.tsv` - mtDNA variant/heteroplasmy calls
  - `*.depthTable.txt` - Per-position mitochondrial coverage

</details>



### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
