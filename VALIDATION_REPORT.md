# atacportray — Lint & Stub Validation Report

Pipeline: **atacportray/atacportray** v1.0.0dev
Template: nf-core tools 4.0.2 (`is_nfcore: false`, standalone)
Validation host: macOS with Docker Desktop for local integration checks

## Component inventory

| Component type | Count |
|---|---|
| Local custom modules | 7 (ROSE, TelomereHunter, mgatk, QDNAseq, TOBIAS ×3, consensus_peaks, maftools_oncoplot) |
| nf-core registry modules | 33 |
| Local subworkflows | 6 (epigenome, variants, CNV, telomere, mito, utils) |
| nf-test files | 48 |

## `nf-core pipelines lint` (pipeline-level tests)

Run with:

```bash
nf-core pipelines lint
```

Current status after the latest cleanup pass:

**Result: 445 PASS · 23 IGNORED · 40 WARN · 0 FAIL**

Remaining warnings are tracked cleanup items:

- `readme`: `zenodo.XXXXXXX` DOI placeholder, expected until the first Zenodo release.
- registry module/subworkflow version warnings, to be handled by a controlled `nf-core modules update` pass.
- local single-module wrapper subworkflows, kept intentionally for branch-level structure and output naming.

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

> On Genotoul, use the institutional `-profile genotoul` with a fresh work
> directory and real ATAC-seq data for production-scale validation.
