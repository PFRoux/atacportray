# atacportray — Lint & Stub Validation Report

Pipeline: **atacportray/atacportray** v1.0.0dev
Template: nf-core tools 3.5.0 (`is_nfcore: false`, standalone)
Validation host: macOS, no container engine (stub mode)

## Component inventory

| Component type | Count |
|---|---|
| Local custom modules | 7 (ROSE, TelomereHunter, mgatk, QDNAseq, TOBIAS ×3, consensus_peaks, maftools_oncoplot) |
| nf-core registry modules | 33 |
| Local subworkflows | 6 (epigenome, variants, CNV, telomere, mito, utils) |
| nf-test files | 48 |

## `nf-core pipelines lint` (pipeline-level tests)

Run via the Python lint API (the CLI's module-lint step requires cloning
`github.com/nf-core/modules.git`, which the validation sandbox cannot do —
git cannot create `.git` directories here). 26 of 27 pipeline-level tests
executed; `actions_schema_validation` needs network to the GitHub schema and
was skipped.

**Result: 193 PASS · 10 WARN · 4 FAIL**

### The 4 FAILs are false positives (one upstream file)

All four are `merge_markers` hits on
`modules/nf-core/bowtie2/align/tests/main.nf.test.snap`. The flagged
`<<<<<<<` strings are **FASTQ base-quality characters inside a SAM record**
(`AEEAA<<<AAAE/E/AA<<<<<<<<...`), not git conflict markers. This is an
upstream nf-core module test-snapshot, unmodified by this pipeline, and the
`merge_markers` heuristic misfires on quality strings. No action required.

### The 10 WARNs are cosmetic template placeholders

- 1 × `readme`: `zenodo.XXXXXXX` DOI placeholder — resolved after first Zenodo release.
- 9 × `pipeline_todos`: standard nf-core template `TODO nf-core:` comments in
  `base.config`, `nextflow.config`, `main.nf`, `tests/nextflow.config` and
  `methods_description_template.yml` (e.g. "Check the defaults for all
  processes", "Specify any additional parameters here"). These are inert
  guidance comments carried by every generated nf-core pipeline.

## Stub run — wiring validation

`nextflow run . -profile <p> -stub` compiles the DAG and executes every
process body's stub block, proving channel wiring end-to-end without real
tools or data.

| Profile | Branches | Result |
|---|---|---|
| `test` | epigenome only | **Pipeline completed successfully** |
| `test_full` | all five branches | All processes submit in correct topological order |

`test` branch process chain (validated):
FASTP → BOWTIE2_BUILD → BOWTIE2_ALIGN → SAMTOOLS_VIEW/SORT/INDEX →
MACS3_CALLPEAK → DEEPTOOLS_BAMCOVERAGE → ROSE → MULTIQC.

`test_full` additionally submits (validated wiring):
- Variants: BWA_MEM → GATK4_MARKDUPLICATES → BASERECALIBRATOR → APPLYBQSR →
  INDEX_ANALYSIS → {HAPLOTYPECALLER, DEEPVARIANT, FREEBAYES, BCFTOOLS_MPILEUP}
  → ENSEMBLVEP_VEP → (VCF2MAF → MAFTOOLS_ONCOPLOT).
- CNV: QDNASEQ · Telomere: TELOMEREHUNTER · Mito: MGATK — all consuming the
  shared analysis-ready BAM emitted by the variant branch.
- Epigenome footprinting: CONSENSUS_PEAKS, TOBIAS_ATACORRECT/SCOREBIGWIG/BINDETECT.

### Bugs found and fixed during stub validation

1. **Blacklist toggle** — when no ENCODE blacklist is supplied the `combine`
   produced a 2-tuple that broke the 3-arg `.map`. Added an `apply_blacklist`
   boolean to the epigenome subworkflow so filtering is conditional.
2. **APPLYBQSR output suffix** — GATK ApplyBQSR defaults to `.cram`, so
   `emit: bam` was empty and severed the entire downstream (callers + CNV +
   telomere + mito). Set `ext.suffix = 'bam'` in `conf/modules.config`.
3. **`modules_testdata_base_path`** missing from schema — added to
   `nextflow_schema.json`.

## nf-test — custom local modules

`nf-test test modules/local/` (stub), JVM pointed at the sandbox HTTP proxy:

**All 10 custom-module tests PASSED**
(consensus_peaks, maftools_oncoplot, mgatk, qdnaseq, rose, telomerehunter,
tobias/atacorrect, tobias/bindetect, tobias/scorebigwig — plus the epigenome
default test).

## Reproducing this validation

```bash
source nfenv.sh                 # sets JAVA_HOME, XDG dirs, JVM proxy, git redirection
cd build/atacportray
# Stub wiring:
nextflow run . -profile test       --input <ss.csv> --fasta <genome.fa> --outdir out -stub
nextflow run . -profile test_full  --input <ss.csv> --fasta <genome.fa> --outdir out -stub
# Local-module unit tests:
nf-test test modules/local/
```

> On a container-capable Linux host with network, run the full
> `nf-core pipelines lint` and `-profile test,docker` (non-stub) for
> real-data validation; the module-lint step and `actions_schema_validation`
> will then execute too.
