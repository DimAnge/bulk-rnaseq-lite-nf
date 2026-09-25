# bulk-rnaseq-lite-nf — Design

Date: 2026-09-25
Status: reviewed, round 1 changes applied

## Purpose

A compact Nextflow (DSL2) pipeline for bulk RNA-seq: raw FASTQ → QC → trimming →
alignment → gene counts → differential expression (DESeq2 + edgeR) → GO enrichment →
one self-contained HTML report.

Built for **lab use** and as a **learning resource for Nextflow**. Two audiences:

1. The author's portfolio — the code should be small, readable, and show sound
   Nextflow and RNA-seq practice.
2. A junior bioinformatician inheriting the lab's analyses — they must be able to run a
   full analysis from the docs alone and re-run the downstream report without touching
   Nextflow.

Success criteria:

- `nextflow run . -profile test,docker` completes in minutes on a 4-core / 15 GB laptop.
- A real human or mouse analysis runs on the lab HPC server (no scheduler, local
  executor) with only a samplesheet and `--genome`.
- The Rmd renders both inside the pipeline and standalone in RStudio.
- The repo is private on GitHub (`bulk-rnaseq-lite-nf`) until finished.

## Tools

| Step | Tool | Notes |
|---|---|---|
| Raw QC | FastQC | on raw reads |
| Trimming | fastp | auto adapter detection, poly-G trimming, JSON for MultiQC |
| Alignment | STAR | `--quantMode GeneCounts` also used for strand detection |
| Counting | featureCounts (Subread) | one matrix for all samples |
| QC aggregation | MultiQC | FastQC, fastp, STAR, featureCounts |
| DE | DESeq2, edgeR | same filtered counts to both |
| Enrichment | clusterProfiler (GO ORA) | offline via `org.Hs.eg.db` / `org.Mm.eg.db` |
| Report | R Markdown | self-contained HTML |

fastp was chosen over Trimmomatic: faster, auto adapter detection (no adapter FASTA),
poly-G handling for NovaSeq data, native MultiQC support.

## Repository layout

```
bulk-rnaseq-lite-nf/
├── main.nf                  # workflow wiring + input/param checks
├── nextflow.config          # defaults, profiles, resources, containers
├── params.yaml              # example user config (-params-file params.yaml)
├── modules/
│   ├── download_refs.nf     # Ensembl FASTA + GTF → storeDir cache
│   ├── star_index.nf        # full or sparse index → storeDir cache
│   ├── fastqc.nf
│   ├── fastp.nf
│   ├── star_align.nf
│   ├── infer_strand.nf      # picks strandedness from ReadsPerGene.out.tab
│   ├── featurecounts.nf
│   ├── multiqc.nf
│   └── report.nf            # renders the Rmd
├── bin/                     # small helper scripts (infer_strand, gene_names from GTF)
├── report/
│   └── rnaseq_report.Rmd
├── docker/
│   └── Dockerfile           # R report image
├── assets/test/             # tiny FASTQs, samplesheet, chr22 slice FASTA + GTF
├── .github/workflows/ci.yml
├── README.md
└── docs/
    ├── usage.md
    └── output.md
```

Plain DSL2, one module per tool. No nf-core template, no nf-schema plugin.

## Data flow

```
samplesheet.csv ─► FASTQC ───────────────────────────────┐
                └► FASTP ─► STAR_ALIGN ─► INFER_STRAND ─┐│
--genome ─► DOWNLOAD_REFS ─► STAR_INDEX ─┘             ▼▼
                         GTF ────────────► FEATURECOUNTS ─► MULTIQC
                                              │
            counts + samplesheet + gene names + QC summary ─► REPORT → report.html
```

- With `--skip_trimming`, STAR takes raw reads.
- `INFER_STRAND` runs per sample; featureCounts runs once with a single strandedness
  value. If samples disagree, the pipeline stops with a message listing per-sample
  calls (mixed library types in one run is almost always an error).
- When `--fasta`/`--gtf`/`--star_index` are supplied, the corresponding download/index
  steps are skipped.

## Input: samplesheet

```
sample,fastq_1,fastq_2,condition,batch
ctrl_1,/data/ctrl_1_R1.fq.gz,/data/ctrl_1_R2.fq.gz,control,A
trt_1,/data/trt_1.fq.gz,,treated,A
```

- Required columns: `sample`, `fastq_1`, `condition`.
- Empty `fastq_2` → single-end for that sample. No global paired/single flag.
- `batch` and any extra columns are optional and usable in `design`.

## Parameters

Set in `params.yaml` (`-params-file`) or on the command line; command line wins.

| Param | Default | Notes |
|---|---|---|
| `input` | required | samplesheet path |
| `genome` | `GRCh38` | `GRCh38` or `GRCm39` |
| `ensembl_release` | `112` | pinned; new value → new cache folder |
| `fasta`, `gtf`, `star_index` | null | overrides |
| `genome_cache` | `~/.rnaseq-refs` | persistent reference cache |
| `strandedness` | `auto` | `auto`, `forward`, `reverse`, `unstranded` |
| `read_length` | `100` | STAR `--sjdbOverhang` = read_length − 1 |
| `design` | `~ condition` | e.g. `~ batch + condition`; last term must be `condition` |
| `contrasts` | null | `"treated_vs_control,KO_vs_WT"`; null → every level vs the first level |
| `padj_cutoff` | `0.05` | both methods |
| `lfc_cutoff` | `0.58` | log2 (≈1.5-fold), both methods |
| `skip_trimming` | false | |
| `skip_enrichment` | false | |
| `outdir` | `results` | |
| `max_cpus` | `8` | resource cap |
| `max_memory` | `32.GB` | resource cap |

## Profiles

- `docker` — containers per process, local executor. Normal use.
- `low_memory` — sparse STAR index (`--genomeSAsparseD 3`, `--limitGenomeGenerateRAM`
  sized to fit ~16 GB). Combined as `-profile docker,low_memory`.
- `test` — bundled test data, small resource caps. `-profile test,docker`.

No SLURM profile in v1; adding one later is a single profile block.

## Reference cache

```
~/.rnaseq-refs/<genome>/ensembl_<release>/
  ├── genome.fa
  ├── genes.gtf
  ├── star_index/           # full
  └── star_index_sparse3/   # low_memory
```

- Download: Ensembl primary assembly FASTA + matching GTF for the pinned release.
- `storeDir` makes both download and index steps no-ops when outputs exist, so each
  runs once per genome/release/index type.
- Docs state the one-time cost for human: ~1 GB download, ~30 GB disk and ~1 h for the
  full index.

## Strandedness auto-detection

STAR `ReadsPerGene.out.tab` columns 2/3/4 = unstranded / forward / reverse counts.

- Fraction = forward / (forward + reverse) using column 3 and 4 totals (gene rows only).
  - ≥ 0.8 → `forward`; ≤ 0.2 → `reverse` (stranded).
  - 0.4–0.6 → `unstranded` (random orientation, ~50/50).
  - Anything else (0.2–0.4 or 0.6–0.8) → ambiguous: pipeline stops and asks for an
    explicit `--strandedness`, printing per-sample fractions.
  - Thresholds are named constants in `bin/infer_strand`.
- Mapped to featureCounts `-s 1 / 2 / 0`.
- Detected value is written to the report and MultiQC.
- Explicit `--strandedness` skips detection.

## Containers

- Tools: pinned BioContainers images (`quay.io/biocontainers/<tool>:<version>--<build>`).
- Report: custom image from `docker/Dockerfile` (base `rocker/r-ver` or
  `bioconductor/bioconductor_docker`, pinned) with rmarkdown, DESeq2, edgeR, limma,
  apeglm, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, ggplot2, pheatmap, DT.
- Built locally: `docker build -t bulk-rnaseq-lite-report docker/`. Publishing to a
  registry is optional, later.

## Report (report/rnaseq_report.Rmd)

Inputs: `counts.tsv`, `samplesheet.csv`, `gene_names.tsv` (Ensembl ID → symbol, parsed
from the GTF), `qc_summary.tsv`, and the `params:` header (design, contrasts, cutoffs,
genome, detected strandedness, skip_enrichment).

Every chunk opens with a plain-English comment explaining what it does and why.

Sections:

1. **Run summary** — parameters, genome + Ensembl release used, detected strandedness,
   tool versions, sample table.
2. **Upstream QC** — per-sample raw reads, % retained after fastp, % uniquely mapped,
   % assigned by featureCounts; low values flagged. Link to MultiQC.
3. **Count QC** — library sizes, genes detected, log-CPM distributions, genes removed by
   `edgeR::filterByExpr`.
4. **Sample relationships** — PCA on VST (colour = condition, shape = batch), PC1–4
   scree, sample-distance heatmap; if `batch` in design, PCA before/after
   `limma::removeBatchEffect` (visual only, not used for DE).
5. **Differential expression**, one tab per contrast — DESeq2 (Wald, apeglm shrinkage)
   and edgeR (quasi-likelihood F-test) side by side: MA plot, volcano, top-50 heatmap,
   searchable `DT` table.
6. **DESeq2 vs edgeR concordance** — overlap of significant genes, log2FC scatter.
7. **GO enrichment** (unless skipped) — clusterProfiler `enrichGO` (BP/MF/CC) on DESeq2
   significant genes only (padj < `padj_cutoff` and |log2FC| > `lfc_cutoff`), up and down
   separately; dot plots + tables. Background = all
   genes after filtering.
8. **Exports** — list of written CSVs.

Outputs: self-contained `report.html` and
`tables/<contrast>_{deseq2,edger,go_up,go_down}.csv`.

Behaviour rules:

- Both methods use the same filtered count matrix.
- Contrast `A_vs_B` → A is numerator. Unknown level → stop with valid levels listed.
- Any group with < 2 replicates → warning in report; still runs.
- Batch fully confounded with condition → batch dropped from design, warning in report.

Standalone: the same Rmd knits in RStudio from a `results/` folder by editing the
`params:` header, or via
`docker run ... Rscript -e "rmarkdown::render('rnaseq_report.Rmd', params=list(...))"`.

## Outputs

```
results/
├── fastqc/  fastp/  star/  featurecounts/  multiqc/
├── report/  report.html, tables/
└── pipeline_info/  timeline.html, trace.txt, report.html
```

## Error handling

- Pre-flight checks in `main.nf` before any process runs: required samplesheet columns,
  FASTQ files exist, unique sample names, valid `genome` and `strandedness`, every
  variable in `design` exists as a samplesheet column.
- STAR and featureCounts: retry once with doubled memory (capped by `max_memory`), then
  stop.
- Low featureCounts assignment rate (< 50%) flagged in report as likely wrong
  strandedness or annotation.
- `-resume` documented prominently.

## Testing

- `-profile test,docker`: 4 subsampled human samples (2 control, 2 treated; mix of
  paired-end and single-end) against a chr22 slice FASTA + GTF in `assets/test/`.
  Expected runtime: minutes on the dev laptop.
- GitHub Actions: runs the test profile on push (includes report render).
- `bin/` helper scripts each carry a small self-check (e.g. strand inference on a
  fixture `ReadsPerGene.out.tab`).

## Documentation

- `README.md` — "Built for lab use and as a learning resource for Nextflow"; what it
  does, flow diagram, 3-command quickstart, citations.
- `docs/usage.md` — all parameters, references and cache, low-memory mode, re-running
  only the report, troubleshooting (low mapping, wrong strand, OOM).
- `docs/output.md` — every output file and how to read each report plot.

## Out of scope for v1

- SLURM or other schedulers.
- KEGG / GSEA / other enrichment databases.
- Species beyond human and mouse (possible via `--fasta`/`--gtf` but no enrichment).
- Singularity/Apptainer profile.
- Pseudo-alignment (Salmon/kallisto).

## Repository

- GitHub: `bulk-rnaseq-lite-nf`, private until finished. Created via `gh` during
  implementation.
