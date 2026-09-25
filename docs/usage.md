# Usage

## Running

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --genome GRCh38
```

Parameters can go on the command line (`--name value`) or in a YAML file
(`-params-file params.yaml`, see the example in the repo root). The command line wins.

**If a run stops** (error, power cut, closed laptop), fix the cause and add `-resume`:
finished steps are reused and only the rest runs again.

```bash
nextflow run main.nf -profile docker -params-file params.yaml -resume
```

## Samplesheet

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Unique name: letters, digits, `_`, `.`, `-` |
| `fastq_1` | yes | Read 1 (or the only read), gzipped: `.fastq.gz` / `.fq.gz` |
| `fastq_2` | no | Read 2; leave empty for single-end |
| `condition` | yes | Group to compare. Starts with a letter; letters, digits, `_`, `.` only; must not contain `_vs_` |
| any other | no | e.g. `batch`, `donor`, `RIN`; usable in `--design` |

Relative FASTQ paths are read relative to the samplesheet's folder. Files saved from Excel are fine.
The first condition that appears is the reference level.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--input` | required | Samplesheet CSV |
| `--outdir` | `results` | Output folder |
| `--genome` | `GRCh38` | `GRCh38` (human) or `GRCm39` (mouse) |
| `--ensembl_release` | `112` | Ensembl release to download |
| `--fasta`, `--gtf` | none | Your own genome + annotation (give both). Any species; GO only for human/mouse |
| `--star_index` | none | Existing STAR index folder (must have been built with a GTF) |
| `--genome_cache` | `~/.rnaseq-refs` | Where downloaded references and indexes are kept |
| `--strandedness` | `auto` | `auto`, `forward`, `reverse` or `unstranded` |
| `--read_length` | `100` | Sets STAR `--sjdbOverhang` when building an index; 100 suits most data |
| `--skip_trimming` | `false` | Align raw reads without fastp |
| `--save_bam` | `false` | Copy BAM files to `results/star/` (they are large) |
| `--design` | `~ condition` | Model formula; must end in `condition`, e.g. `~ batch + condition` |
| `--contrasts` | all vs reference | Comma-separated `A_vs_B` list, e.g. `treated_vs_control,KO_vs_WT` |
| `--padj_cutoff` | `0.05` | Adjusted p-value cutoff (DESeq2 and edgeR) |
| `--lfc_cutoff` | `0.58` | Absolute log2 fold-change cutoff (0.58 ≈ 1.5-fold) |
| `--skip_enrichment` | `false` | Skip GO enrichment |
| `--max_cpus` | `8` | No task uses more CPUs than this |
| `--max_memory` | `32.GB` | No task uses more memory than this |

## Profiles

| Profile | Use |
|---|---|
| `docker` | Normal runs (always include it) |
| `low_memory` | Machines with ~16 GB RAM: sparse STAR index, slower alignment. `-profile docker,low_memory` |
| `test` | Bundled tiny dataset: `-profile test,docker` |

## Reference genomes and the cache

With `--genome GRCh38` or `GRCm39`, the first run downloads the Ensembl primary-assembly FASTA
and GTF for `--ensembl_release`, then builds the STAR index, all under:

```
~/.rnaseq-refs/<genome>/ensembl_<release>/
  genome.fa  genes.gtf  star_index/  star_index_sparse3/   (the last one only with low_memory)
```

Nextflow checks this folder first and skips the download/index steps if the files are there,
so this happens **once per genome and release**. One-time cost for human: ~1 GB download,
~30 GB disk, ~1 h and ~32 GB RAM for the full index (~16 GB for the sparse one).

To share one cache between users on a server, point everyone at the same folder:
`--genome_cache /shared/refs/rnaseq`.

With `--fasta` and `--gtf` the index is built into `<outdir>/reference/`. Reuse it in later runs
with `--star_index <outdir>/reference/star_index --gtf <your.gtf>`.

## Strandedness

`auto` looks at STAR's per-gene counts for each sample: the fraction of reads on the forward strand.

| Forward fraction | Call |
|---|---|
| ≥ 0.8 | `forward` |
| ≤ 0.2 | `reverse` (most dUTP/TruSeq stranded kits) |
| 0.4–0.6 | `unstranded` |
| anything else | the run stops and asks you to set `--strandedness` |

The run also stops if samples disagree. The per-sample numbers are in
`results/star/strandedness_mqc.tsv`, MultiQC and the report.

## Re-running only the report

Changing a contrast, cutoff or design doesn't require re-running alignment.

**Option A: through the pipeline.** Run the same command with `-resume` and the new parameter;
only the report step reruns.

```bash
nextflow run main.nf -profile docker -params-file params.yaml -resume --contrasts KO_vs_WT
```

**Option B: RStudio.** Open `report/rnaseq_report.Rmd`. The `params:` block at the top points at
`../results/...`, which works when your `results/` folder is next to `report/`; otherwise edit the paths.
Change what you need and click **Knit**. Needs R with the packages listed in `docker/Dockerfile`.

**Option C: Docker, no local R needed.** From the repo folder:

```bash
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/w -w /w/report bulk-rnaseq-lite-report:1.0 \
  Rscript -e 'rmarkdown::render("rnaseq_report.Rmd", params = list(contrasts = "KO_vs_WT", lfc_cutoff = 1))'
```

The report is written to `report/rnaseq_report.html` and tables to `report/tables/`.

## Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `Strandedness is ambiguous` / `Samples disagree` | Mixed library types or odd data. Check `results/star/strandedness_mqc.tsv`; set `--strandedness` if you know the kit. |
| Low "% assigned to genes" warning | Wrong strandedness (try another `--strandedness`), or a GTF that doesn't match the genome. |
| Low "% uniquely mapped" | Wrong species, contamination or rRNA; look at FastQC/MultiQC. |
| STAR killed / exit 137 | Not enough RAM. Use `-profile docker,low_memory` or a bigger machine; lower `--max_cpus` if many STAR jobs run at once. |
| `--quantMode GeneCounts` error with `--star_index` | The index was built without a GTF. Rebuild it (drop `--star_index`). |
| Docker permission denied | Add your user to the `docker` group or ask your admin. |
| `Unable to find image 'bulk-rnaseq-lite-report:1.0'` | Build it: `docker build -t bulk-rnaseq-lite-report:1.0 docker/`. |
