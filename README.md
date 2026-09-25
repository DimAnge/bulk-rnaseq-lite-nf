# bulk-rnaseq-lite-nf

A compact Nextflow pipeline for bulk RNA-seq: raw FASTQ files in, gene counts,
differential expression (DESeq2 **and** edgeR), GO enrichment and one HTML report out.

Built for lab use and as a learning resource for Nextflow: every step is one short,
commented file in `modules/`, and `main.nf` reads top to bottom.

## What it does

```
samplesheet.csv ─► FastQC
                └► fastp ─► STAR ─► strandedness check ─► featureCounts ─► MultiQC
                                                             │
                                  DESeq2 + edgeR + GO ◄──────┘ ─► report.html
```

- Human (GRCh38) and mouse (GRCm39), references downloaded from Ensembl **once** and cached.
- Single- and paired-end samples, even mixed in one run.
- Strandedness detected automatically (or set it yourself).
- Design with batch or other covariates: `--design '~ batch + condition'`.
- Report: QC tables, PCA, sample distances, MA/volcano plots, heatmaps, DESeq2 vs edgeR comparison, GO dot plots, CSV tables.

## Quickstart

Requirements: [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 25.04 and Docker.

```bash
# 1. Build the R report image (once, 10-30 min)
docker build -t bulk-rnaseq-lite-report:1.0 docker/

# 2. Run the bundled test data (a few minutes)
nextflow run main.nf -profile test,docker
#    -> open results/report/report.html

# 3. Run your own data
nextflow run main.nf -profile docker --input samplesheet.csv --genome GRCh38
```

The first real run for a genome downloads it (~1 GB) and builds the STAR index
(~30 GB disk, ~1 h, ~32 GB RAM) into `~/.rnaseq-refs/`. Later runs reuse it.
On a machine with ~16 GB RAM add the low-memory profile: `-profile docker,low_memory`.

## Samplesheet

```csv
sample,fastq_1,fastq_2,condition,batch
ctrl_1,/data/ctrl_1_R1.fastq.gz,/data/ctrl_1_R2.fastq.gz,control,A
trt_1,/data/trt_1_R1.fastq.gz,,treated,A
```

Leave `fastq_2` empty for single-end. `batch` (and any other column) is optional.
The first condition listed is the reference.

## Documentation

- [docs/usage.md](docs/usage.md): all parameters, references and cache, low-memory mode, re-running only the report, troubleshooting.
- [docs/output.md](docs/output.md): every output file and how to read the report.

## Tools and citations

FastQC · fastp (Chen et al. 2018, *Bioinformatics*) · STAR (Dobin et al. 2013, *Bioinformatics*) ·
featureCounts (Liao et al. 2014, *Bioinformatics*) · MultiQC (Ewels et al. 2016, *Bioinformatics*) ·
DESeq2 (Love et al. 2014, *Genome Biology*) · apeglm (Zhu et al. 2019, *Bioinformatics*) ·
edgeR (Chen et al. 2025, *Nucleic Acids Research*) · clusterProfiler (Xu et al. 2024, *Nature Protocols*) ·
Nextflow (Di Tommaso et al. 2017, *Nature Biotechnology*).
