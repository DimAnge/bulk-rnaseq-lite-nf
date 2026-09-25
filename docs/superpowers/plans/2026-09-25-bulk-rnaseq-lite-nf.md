# bulk-rnaseq-lite-nf Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A compact Nextflow DSL2 pipeline that takes bulk RNA-seq FASTQs to counts, DESeq2 + edgeR differential expression, GO enrichment and one self-contained HTML report, running in Docker for human and mouse.

**Architecture:** One process per tool under `modules/`, wired in `main.nf`, which also validates all inputs before any process starts. Small Python helpers in `bin/` handle strand inference, count merging, gene names and QC summaries. A single R Markdown report (plus one child template) runs both inside the pipeline and standalone in RStudio. References download once into a persistent `storeDir` cache.

**Tech Stack:** Nextflow ≥ 25.04 (strict-syntax compatible), Docker, FastQC 0.12.1, fastp 0.23.4, STAR 2.7.11b, Subread/featureCounts 2.0.6, MultiQC 1.25.1, Python 3.12 (stdlib only), R/Bioconductor 3.22 (DESeq2, edgeR, limma, apeglm, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, rmarkdown, ggplot2, pheatmap, DT).

**Spec:** `docs/superpowers/specs/2026-09-25-bulk-rnaseq-lite-nf-design.md`

## Global Constraints

- Repo name `bulk-rnaseq-lite-nf`; GitHub repo **private** until finished.
- README framing: "Built for lab use and as a learning resource for Nextflow". No nf-core comparison.
- Plain DSL2, one module per tool; no nf-core template, no nf-schema plugin.
- Code must pass `nextflow lint` (strict syntax): no `for`/`while` loops in Nextflow code, explicit closure parameters (no implicit `it`), `error()` instead of `exit`, lowercase `channel`.
- Genomes: `GRCh38` (human) and `GRCm39` (mouse), Ensembl release default `112`.
- Default cutoffs: `padj_cutoff = 0.05`, `lfc_cutoff = 0.58`; used by both DESeq2 and edgeR.
- Strandedness default `auto`; bands: ≥ 0.8 forward, ≤ 0.2 reverse, 0.4–0.6 unstranded, otherwise stop and ask.
- Reference cache default `~/.rnaseq-refs/<genome>/ensembl_<release>/`; never re-downloads or re-indexes if present.
- Report shows genome, Ensembl release and detected strandedness. GO runs only on significant DESeq2 genes (up and down separately).
- Every R chunk starts with a plain-English comment for a junior reader.
- Python helpers: stdlib only.
- Commits end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. **Samplesheet saved from Excel (UTF-8 BOM + CRLF line endings)** → parses normally, no "missing column sample" error. Test: Task 4 validation case `bom_crlf`.
2. **Every sample's FASTQs have the same basename** (`reads/<sample>/R1.fastq.gz`) → no staging or output collisions. Test: Task 3 test data uses exactly this layout; Task 5 run must pass.
3. **Nothing significant / fewer than 1000 genes** → report still renders with explanatory messages (VST fallback, no heatmap, no GO). Test: Task 7 standalone render with `lfc_cutoff = 100`.
4. **Batch column confounded with condition** → design falls back to `~ condition` with a visible warning, report renders. Test: Task 7 standalone render with batch == condition.
5. **Second run with the same reference** → STAR index is reused, not rebuilt. Test: Task 5 re-run and index mtime unchanged.

---

## File map

```
bulk-rnaseq-lite-nf/
├── .gitignore
├── .github/workflows/ci.yml          # Task 9
├── README.md                          # Task 8
├── main.nf                            # Task 4 (checks), Task 5 (upstream wiring), Task 7 (report wiring)
├── nextflow.config                    # Task 4
├── params.yaml                        # Task 8
├── assets/NO_FILE                     # Task 5 (placeholder input)
├── assets/test/                       # Task 3 (generator + generated data)
├── bin/infer_strand.py                # Task 2
├── bin/gene_names.py                  # Task 2
├── bin/merge_counts.py                # Task 2
├── bin/qc_summary.py                  # Task 2
├── modules/download_refs.nf           # Task 5
├── modules/star_index.nf              # Task 5
├── modules/fastqc.nf                  # Task 5
├── modules/fastp.nf                   # Task 5
├── modules/star_align.nf              # Task 5
├── modules/infer_strand.nf            # Task 5
├── modules/featurecounts.nf           # Task 5
├── modules/helpers.nf                 # Task 5 (MERGE_COUNTS, GENE_NAMES, QC_SUMMARY)
├── modules/multiqc.nf                 # Task 5
├── modules/report.nf                  # Task 7
├── docker/Dockerfile                  # Task 6
├── report/rnaseq_report.Rmd           # Task 7
├── report/_contrast.Rmd               # Task 7 (per-contrast child template)
├── tests/test_bin.py                  # Task 2
├── tests/validation.sh                # Task 4
├── tests/check_results.py             # Task 7
└── docs/usage.md, docs/output.md      # Task 8
```

Two small deviations from the spec's file list, both for simplicity: featureCounts runs **per sample** (a single featureCounts call cannot mix paired- and single-end BAMs) and `MERGE_COUNTS` builds the single matrix; the three tiny Python-wrapper processes share `modules/helpers.nf`.

---

### Task 1: Repository skeleton and private GitHub repo

**Files:**
- Create: `.gitignore`

**Interfaces:**
- Consumes: existing local repo at `/home/ange/projects/bulk-rnaseq-lite-nf` (branch `main`, spec committed).
- Produces: GitHub remote `origin` → private `DimAnge/bulk-rnaseq-lite-nf`.

- [ ] **Step 1: Write `.gitignore`**

```gitignore
# Nextflow
work/
.nextflow/
.nextflow.log*
results/
# Test-data downloads (generator cache)
assets/test/_downloads/
# R / RStudio
.Rproj.user/
.Rhistory
report/*.html
report/tables/
report/*_files/
report/*.knit.md
# Python
__pycache__/
```

- [ ] **Step 2: Commit**

```bash
cd /home/ange/projects/bulk-rnaseq-lite-nf
git add .gitignore
git commit -m "chore: add .gitignore

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 3: Create the private GitHub repo and push**

Run: `gh repo create bulk-rnaseq-lite-nf --private --source . --remote origin --push --description "Compact Nextflow bulk RNA-seq pipeline: FASTQ to DESeq2/edgeR report"`
Expected: prints the repo URL; `gh repo view --json visibility -q .visibility` prints `PRIVATE`.

---

### Task 2: Python helper scripts

**Files:**
- Create: `bin/infer_strand.py`, `bin/gene_names.py`, `bin/merge_counts.py`, `bin/qc_summary.py`
- Test: `tests/test_bin.py`

**Interfaces:**
- Produces (CLI, all write into the current directory):
  - `infer_strand.py <id>.ReadsPerGene.out.tab ...` → `strandedness.txt` (one word: `forward|reverse|unstranded`), `strandedness_mqc.tsv` (columns `Sample, forward_fraction, call`, MultiQC comment header). Exits non-zero with a per-sample table if calls are ambiguous or disagree.
  - `gene_names.py <genes.gtf[.gz]>` → `gene_names.tsv` (columns `gene_id, gene_name, gene_biotype`).
  - `merge_counts.py <id>.featureCounts.txt ...` → `counts.tsv` (columns `gene_id, <sample>...`, samples sorted by name).
  - `qc_summary.py <files...>` accepting `<id>.fastp.json`, `<id>.Log.final.out`, `<id>.featureCounts.txt.summary` → `qc_summary.tsv` (columns `sample, raw_reads, trimmed_reads, pct_retained, star_input_reads, pct_uniquely_mapped, pct_assigned`; missing values `NA`).
  - Python functions used by tests: `infer_strand.call(frac) -> str`, `infer_strand.main(paths)`, `gene_names.gene_table(path) -> dict`, `gene_names.main(path)`, `merge_counts.main(paths)`, `qc_summary.main(paths)`.

- [ ] **Step 1: Write the failing tests**

`tests/test_bin.py`:

```python
#!/usr/bin/env python3
"""Self-checks for the helper scripts in bin/. Run: python3 tests/test_bin.py"""
import gzip
import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
import gene_names  # noqa: E402
import infer_strand  # noqa: E402
import merge_counts  # noqa: E402
import qc_summary  # noqa: E402


def write(path, text):
    Path(path).write_text(text)
    return path


def reads_per_gene(fwd, rev):
    return ("N_unmapped\t10\t10\t10\nN_multimapping\t5\t5\t5\n"
            "N_noFeature\t3\t40\t40\nN_ambiguous\t1\t1\t1\n"
            f"GENE1\t{fwd + rev}\t{fwd}\t{rev}\nGENE2\t0\t0\t0\n")


def expect_exit(fn, text):
    try:
        fn()
    except SystemExit as e:
        assert text in str(e), f"expected '{text}' in: {e}"
        return
    raise AssertionError(f"expected SystemExit containing '{text}'")


def test_infer_strand_call():
    assert infer_strand.call(0.95) == "forward"
    assert infer_strand.call(0.05) == "reverse"
    assert infer_strand.call(0.20) == "reverse"
    assert infer_strand.call(0.50) == "unstranded"
    assert infer_strand.call(0.70) == "ambiguous"
    assert infer_strand.call(0.30) == "ambiguous"


def test_infer_strand_main_agreeing_samples():
    a = write("a.ReadsPerGene.out.tab", reads_per_gene(5, 95))
    b = write("b.ReadsPerGene.out.tab", reads_per_gene(3, 97))
    infer_strand.main([a, b])
    assert Path("strandedness.txt").read_text().strip() == "reverse"
    mqc = Path("strandedness_mqc.tsv").read_text()
    assert "# plot_type: 'table'" in mqc
    assert "a\t0.050\treverse" in mqc


def test_infer_strand_main_disagreeing_samples():
    a = write("a.ReadsPerGene.out.tab", reads_per_gene(5, 95))
    c = write("c.ReadsPerGene.out.tab", reads_per_gene(95, 5))
    expect_exit(lambda: infer_strand.main([a, c]), "disagree")
    assert Path("strandedness_mqc.tsv").exists()  # table still written for debugging


def test_infer_strand_main_ambiguous():
    d = write("d.ReadsPerGene.out.tab", reads_per_gene(70, 30))
    expect_exit(lambda: infer_strand.main([d]), "ambiguous")


def test_infer_strand_no_assigned_reads():
    e = write("e.ReadsPerGene.out.tab", reads_per_gene(0, 0))
    expect_exit(lambda: infer_strand.main([e]), "no reads were assigned")


def test_gene_names():
    gtf = ('#!genome-build GRCh38\n'
           '22\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "ABC"; gene_biotype "protein_coding";\n'
           '22\tensembl\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ABC";\n'
           '22\tgencode\texon\t5\t50\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_type "lncRNA";\n')
    with gzip.open("g.gtf.gz", "wt") as fh:
        fh.write(gtf)
    table = gene_names.gene_table("g.gtf.gz")
    assert table == {"G1": ("ABC", "protein_coding"), "G2": ("G2", "lncRNA")}
    gene_names.main("g.gtf.gz")
    lines = Path("gene_names.tsv").read_text().splitlines()
    assert lines == ["gene_id\tgene_name\tgene_biotype", "G1\tABC\tprotein_coding", "G2\tG2\tlncRNA"]


def fc_file(name, counts):
    body = "".join(f"{g}\t22\t1\t10\t+\t10\t{n}\n" for g, n in counts)
    return write(name, f"# Program:featureCounts v2.0.6\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tx.bam\n{body}")


def test_merge_counts():
    b = fc_file("s2.featureCounts.txt", [("G1", 7), ("G2", 0)])
    a = fc_file("s1.featureCounts.txt", [("G1", 5), ("G2", 3)])
    merge_counts.main([b, a])
    assert Path("counts.tsv").read_text().splitlines() == ["gene_id\ts1\ts2", "G1\t5\t7", "G2\t3\t0"]


def test_merge_counts_mismatched_genes():
    a = fc_file("s1.featureCounts.txt", [("G1", 5)])
    b = fc_file("s2.featureCounts.txt", [("G9", 5)])
    expect_exit(lambda: merge_counts.main([a, b]), "gene list differs")


def test_qc_summary():
    json.dump({"summary": {"before_filtering": {"total_reads": 1000},
                           "after_filtering": {"total_reads": 900}}}, open("s1.fastp.json", "w"))
    star_log = ("                          Number of input reads |\t500\n"
                "                        Uniquely mapped reads % |\t92.50%\n")
    write("s1.Log.final.out", star_log)
    write("s2.Log.final.out", star_log)
    fc = "Status\ts1.bam\nAssigned\t300\nUnassigned_NoFeatures\t100\nUnassigned_Unmapped\t0\n"
    write("s1.featureCounts.txt.summary", fc)
    write("s2.featureCounts.txt.summary", fc)
    qc_summary.main(["s1.fastp.json", "s1.Log.final.out", "s1.featureCounts.txt.summary",
                     "s2.Log.final.out", "s2.featureCounts.txt.summary"])
    lines = Path("qc_summary.tsv").read_text().splitlines()
    assert lines[0] == "sample\traw_reads\ttrimmed_reads\tpct_retained\tstar_input_reads\tpct_uniquely_mapped\tpct_assigned"
    assert lines[1] == "s1\t1000\t900\t90.0\t500\t92.5\t75.0"
    assert lines[2] == "s2\tNA\tNA\tNA\t500\t92.5\t75.0"  # --skip_trimming: no fastp JSON


if __name__ == "__main__":
    tests = [(name, fn) for name, fn in sorted(globals().items()) if name.startswith("test_")]
    for name, fn in tests:
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            try:
                fn()
            finally:
                os.chdir(cwd)
        print("ok", name)
    print(f"{len(tests)} tests passed")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 tests/test_bin.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'gene_names'`.

- [ ] **Step 3: Write `bin/infer_strand.py`**

```python
#!/usr/bin/env python3
"""Infer library strandedness from STAR ReadsPerGene.out.tab files.

STAR --quantMode GeneCounts writes, per gene: unstranded (col 2), forward (col 3) and
reverse (col 4) counts. For each sample: fraction = forward / (forward + reverse).
  >= 0.8 forward | <= 0.2 reverse | 0.4-0.6 unstranded | anything else: ambiguous
All samples must agree, otherwise the run stops and asks for --strandedness.

Usage:  infer_strand.py SAMPLE.ReadsPerGene.out.tab [...]
Writes: strandedness.txt (the call), strandedness_mqc.tsv (per-sample table for MultiQC + report)
"""
import os
import sys

FORWARD_MIN = 0.8
REVERSE_MAX = 0.2
UNSTRANDED_LOW, UNSTRANDED_HIGH = 0.4, 0.6
SUFFIX = ".ReadsPerGene.out.tab"

MQC_HEADER = """# id: 'strandedness'
# section_name: 'Strandedness'
# description: 'Fraction of gene-assigned reads on the forward strand (STAR GeneCounts). ~0.5 = unstranded, high = forward, low = reverse.'
# plot_type: 'table'
"""


def forward_fraction(path):
    fwd = rev = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("N_"):  # summary rows, not genes
                continue
            cols = line.rstrip("\n").split("\t")
            fwd += int(cols[2])
            rev += int(cols[3])
    if fwd + rev == 0:
        raise SystemExit(f"{path}: no reads were assigned to genes, cannot infer strandedness. "
                         "Check the annotation or set --strandedness.")
    return fwd / (fwd + rev)


def call(frac):
    if frac >= FORWARD_MIN:
        return "forward"
    if frac <= REVERSE_MAX:
        return "reverse"
    if UNSTRANDED_LOW <= frac <= UNSTRANDED_HIGH:
        return "unstranded"
    return "ambiguous"


def main(paths):
    rows = []
    for p in sorted(paths):
        frac = forward_fraction(p)
        rows.append((os.path.basename(p).replace(SUFFIX, ""), frac, call(frac)))
    with open("strandedness_mqc.tsv", "w") as out:
        out.write(MQC_HEADER)
        out.write("Sample\tforward_fraction\tcall\n")
        for sample, frac, c in rows:
            out.write(f"{sample}\t{frac:.3f}\t{c}\n")

    calls = {c for _, _, c in rows}
    table = "\n".join(f"  {s}: forward fraction {f:.3f} -> {c}" for s, f, c in rows)
    if "ambiguous" in calls:
        raise SystemExit("Strandedness is ambiguous for some samples (forward fraction between "
                         "0.2-0.4 or 0.6-0.8). Set --strandedness explicitly.\n" + table)
    if len(calls) > 1:
        raise SystemExit("Samples disagree on strandedness (mixed library types?). "
                         "Set --strandedness explicitly or run them separately.\n" + table)
    with open("strandedness.txt", "w") as out:
        out.write(calls.pop() + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
```

- [ ] **Step 4: Write `bin/gene_names.py`**

```python
#!/usr/bin/env python3
"""Make a gene_id -> gene_name / biotype table from a GTF (plain or .gz).

Works with Ensembl (gene_biotype) and GENCODE (gene_type) GTFs, with or without 'gene'
lines. Genes without a gene_name keep their ID as the name.

Usage:  gene_names.py genes.gtf   ->  gene_names.tsv
"""
import gzip
import re
import sys

ATTR = re.compile(r'(\S+) "([^"]*)"')


def gene_table(path):
    opener = gzip.open if path.endswith(".gz") else open
    genes = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = dict(ATTR.findall(fields[8]))
            gid = attrs.get("gene_id")
            if gid and gid not in genes:
                genes[gid] = (attrs.get("gene_name", gid),
                              attrs.get("gene_biotype", attrs.get("gene_type", "NA")))
    return genes


def main(path, out="gene_names.tsv"):
    with open(out, "w") as fh:
        fh.write("gene_id\tgene_name\tgene_biotype\n")
        for gid, (name, biotype) in gene_table(path).items():
            fh.write(f"{gid}\t{name}\t{biotype}\n")


if __name__ == "__main__":
    main(sys.argv[1])
```

- [ ] **Step 5: Write `bin/merge_counts.py`**

```python
#!/usr/bin/env python3
"""Merge per-sample featureCounts outputs into one gene x sample count matrix.

Usage:  merge_counts.py SAMPLE.featureCounts.txt [...]   ->  counts.tsv
The sample name is the file name without '.featureCounts.txt'.
"""
import os
import sys

SUFFIX = ".featureCounts.txt"


def read_counts(path):
    genes, counts = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid\t"):
                continue
            cols = line.rstrip("\n").split("\t")
            genes.append(cols[0])
            counts.append(cols[6])
    return genes, counts


def main(paths, out="counts.tsv"):
    samples, columns, genes = [], [], None
    for p in sorted(paths, key=os.path.basename):
        g, c = read_counts(p)
        if genes is None:
            genes = g
        elif g != genes:
            raise SystemExit(f"{p}: gene list differs from the other samples (different GTF?)")
        samples.append(os.path.basename(p)[: -len(SUFFIX)])
        columns.append(c)
    with open(out, "w") as fh:
        fh.write("gene_id\t" + "\t".join(samples) + "\n")
        for i, gene in enumerate(genes):
            fh.write(gene + "\t" + "\t".join(col[i] for col in columns) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
```

- [ ] **Step 6: Write `bin/qc_summary.py`**

```python
#!/usr/bin/env python3
"""Collect per-sample QC numbers into one table for the report.

Usage:  qc_summary.py FILES...   ->  qc_summary.tsv
Recognised files (sample name = file name minus suffix):
  SAMPLE.fastp.json                    reads before/after trimming (all mates)
  SAMPLE.Log.final.out                 STAR input reads, % uniquely mapped
  SAMPLE.featureCounts.txt.summary     % of alignments assigned to genes
Missing values (e.g. no fastp with --skip_trimming) are written as NA.
"""
import json
import os
import sys

COLUMNS = ["sample", "raw_reads", "trimmed_reads", "pct_retained",
           "star_input_reads", "pct_uniquely_mapped", "pct_assigned"]


def parse_fastp(path):
    with open(path) as fh:
        s = json.load(fh)["summary"]
    before = s["before_filtering"]["total_reads"]
    after = s["after_filtering"]["total_reads"]
    return {"raw_reads": before, "trimmed_reads": after,
            "pct_retained": round(100 * after / before, 2) if before else "NA"}


def parse_star_log(path):
    vals = {}
    with open(path) as fh:
        for line in fh:
            if "|" in line:
                key, value = line.split("|", 1)
                vals[key.strip()] = value.strip()
    return {"star_input_reads": int(vals["Number of input reads"]),
            "pct_uniquely_mapped": float(vals["Uniquely mapped reads %"].rstrip("%"))}


def parse_fc_summary(path):
    counts = {}
    with open(path) as fh:
        next(fh)  # header: Status <bam>
        for line in fh:
            status, n = line.rstrip("\n").split("\t")
            counts[status] = int(n)
    total = sum(counts.values())
    return {"pct_assigned": round(100 * counts.get("Assigned", 0) / total, 2) if total else "NA"}


PARSERS = {".fastp.json": parse_fastp,
           ".Log.final.out": parse_star_log,
           ".featureCounts.txt.summary": parse_fc_summary}


def main(paths, out="qc_summary.tsv"):
    rows = {}
    for p in paths:
        name = os.path.basename(p)
        for suffix, parser in PARSERS.items():
            if name.endswith(suffix):
                rows.setdefault(name[: -len(suffix)], {}).update(parser(p))
                break
    with open(out, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for sample in sorted(rows):
            values = [sample] + [str(rows[sample].get(c, "NA")) for c in COLUMNS[1:]]
            fh.write("\t".join(values) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
```

- [ ] **Step 7: Make scripts executable and run tests**

Run: `chmod +x bin/*.py && python3 tests/test_bin.py`
Expected: nine `ok test_...` lines, then `9 tests passed`.

- [ ] **Step 8: Commit**

```bash
git add bin/ tests/test_bin.py
git commit -m "feat: add Python helpers for strand inference, counts, gene names, QC summary

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Tiny simulated test dataset

**Files:**
- Create: `assets/test/make_test_data.py`
- Generated (committed): `assets/test/genome.fa`, `assets/test/genes.gtf`, `assets/test/samplesheet.csv`, `assets/test/de_truth.tsv`, `assets/test/reads/<sample>/R1.fastq.gz` (+ `R2.fastq.gz` for paired samples)

**Interfaces:**
- Produces: a 2 Mb slice of chr22 (contig name `22`) with its GTF; 6 samples (3 control, 3 treated; `ctrl_2`/`trt_2` single-end, the rest paired-end); columns `sample,fastq_1,fastq_2,condition,batch` with paths relative to the samplesheet; reverse-stranded (dUTP-like) 50 bp reads; `de_truth.tsv` (`gene_id, direction` with direction `up|down`, 4-fold change).

- [ ] **Step 1: Write `assets/test/make_test_data.py`**

```python
#!/usr/bin/env python3
"""Build the tiny test dataset in assets/test/ (run once; the outputs are committed).

- Downloads chr22 FASTA + the GRCh38 GTF from Ensembl 112 (cached in _downloads/).
- Keeps a 2 Mb slice (chr22:20,000,001-22,000,000) and the genes fully inside it.
- Simulates reverse-stranded (dUTP-like) 50 bp reads for 6 samples. A set of genes is
  4x up and another 4x down in 'treated'; the truth is written to de_truth.tsv.
- Every sample's reads are called R1/R2.fastq.gz in their own folder on purpose: it
  checks the pipeline never relies on FASTQ file names being unique.

Usage: python3 assets/test/make_test_data.py
"""
import gzip
import random
import re
import urllib.request
from pathlib import Path

RELEASE = 112
BASE = f"https://ftp.ensembl.org/pub/release-{RELEASE}"
FASTA_URL = f"{BASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
GTF_URL = f"{BASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{RELEASE}.gtf.gz"
START, END = 20_000_000, 22_000_000  # slice = chr22 positions START+1 .. END
READ_LEN, FRAG_LEN, READS_PER_SAMPLE = 50, 200, 15_000
FOLD, NOISE, MAX_DE = 4.0, 0.15, 20
SEED = 42
SAMPLES = [  # name, condition, batch, paired-end
    ("ctrl_1", "control", "A", True), ("ctrl_2", "control", "B", False), ("ctrl_3", "control", "B", True),
    ("trt_1", "treated", "A", True), ("trt_2", "treated", "B", False), ("trt_3", "treated", "A", True),
]
OUT = Path(__file__).resolve().parent
ATTR = re.compile(r'(\S+) "([^"]*)"')
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def fetch(url):
    dest = OUT / "_downloads" / url.rsplit("/", 1)[1]
    dest.parent.mkdir(exist_ok=True)
    if not dest.exists():
        print("downloading", url)
        urllib.request.urlretrieve(url, dest)
    return dest


def load_slice_sequence(fasta_gz):
    with gzip.open(fasta_gz, "rt") as fh:
        next(fh)  # header
        seq = "".join(line.strip() for line in fh).upper()
    return seq[START:END]


def load_slice_gtf(gtf_gz):
    """GTF records (coordinates shifted to the slice) for chr22 genes fully inside it."""
    keep, records = set(), []
    with gzip.open(gtf_gz, "rt") as fh:
        for line in fh:
            if not line.startswith("22\t"):
                continue
            f = line.rstrip("\n").split("\t")
            attrs = dict(ATTR.findall(f[8]))
            if f[2] == "gene" and int(f[3]) > START and int(f[4]) <= END:
                keep.add(attrs["gene_id"])
            records.append((f, attrs))
    out = []
    for f, attrs in records:
        if attrs["gene_id"] in keep:
            f = f.copy()
            f[3], f[4] = str(int(f[3]) - START), str(int(f[4]) - START)
            out.append((f, attrs))
    return out


def longest_transcripts(records, genome):
    """gene_id -> spliced sequence (sense strand) of its longest protein-coding transcript."""
    exons = {}
    for f, a in records:
        if f[2] == "exon" and a.get("gene_biotype") == "protein_coding":
            exons.setdefault((a["gene_id"], a["transcript_id"], f[6]), []).append((int(f[3]), int(f[4])))
    best = {}
    for (gene, _tx, strand), ex in exons.items():
        seq = "".join(genome[s - 1:e] for s, e in sorted(ex))
        if strand == "-":
            seq = revcomp(seq)
        if len(seq) >= FRAG_LEN and len(seq) > len(best.get(gene, "")):
            best[gene] = seq
    return best


def fastq_gz(path, reads):
    text = "".join(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in reads)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(text.encode(), mtime=0))  # mtime=0: reproducible bytes


def main():
    rng = random.Random(SEED)
    genome = load_slice_sequence(fetch(FASTA_URL))
    records = load_slice_gtf(fetch(GTF_URL))

    (OUT / "genome.fa").write_text(">22\n" + "\n".join(genome[i:i + 60] for i in range(0, len(genome), 60)) + "\n")
    with open(OUT / "genes.gtf", "w") as fh:
        fh.write(f"#!Slice of Ensembl {RELEASE} GRCh38 chr22:{START + 1}-{END}, coordinates shifted\n")
        fh.writelines("\t".join(f) + "\n" for f, _ in records)

    tx = longest_transcripts(records, genome)
    genes = sorted(tx)
    base = {g: rng.lognormvariate(0, 0.5) * len(tx[g]) for g in genes}
    # DE genes come from the better-expressed half so the test has power to find them.
    candidates = sorted(genes, key=lambda g: base[g], reverse=True)[: len(genes) // 2]
    rng.shuffle(candidates)
    n_de = min(MAX_DE, len(candidates) // 2)
    up, down = set(candidates[:n_de]), set(candidates[n_de:2 * n_de])
    with open(OUT / "de_truth.tsv", "w") as fh:
        fh.write("gene_id\tdirection\n")
        fh.writelines(f"{g}\tup\n" for g in sorted(up))
        fh.writelines(f"{g}\tdown\n" for g in sorted(down))

    sheet = ["sample,fastq_1,fastq_2,condition,batch"]
    for name, condition, batch, paired in SAMPLES:
        treated = condition == "treated"
        weights = []
        for g in genes:
            fold = FOLD if (treated and g in up) else (1 / FOLD if (treated and g in down) else 1.0)
            weights.append(base[g] * fold * rng.lognormvariate(0, NOISE))
        r1, r2 = [], []
        for i, g in enumerate(rng.choices(genes, weights=weights, k=READS_PER_SAMPLE)):
            s = rng.randint(0, len(tx[g]) - FRAG_LEN)
            frag = tx[g][s:s + FRAG_LEN]
            r1.append((f"{name}_{i}/1", revcomp(frag)[:READ_LEN]))  # dUTP: read 1 is antisense
            r2.append((f"{name}_{i}/2", frag[:READ_LEN]))
        fastq_gz(OUT / "reads" / name / "R1.fastq.gz", r1)
        if paired:
            fastq_gz(OUT / "reads" / name / "R2.fastq.gz", r2)
        fq2 = f"reads/{name}/R2.fastq.gz" if paired else ""
        sheet.append(f"{name},reads/{name}/R1.fastq.gz,{fq2},{condition},{batch}")
    (OUT / "samplesheet.csv").write_text("\n".join(sheet) + "\n")
    print(f"{len(genes)} simulated genes, {len(up)} up, {len(down)} down")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the generator**

Run: `python3 assets/test/make_test_data.py`
Expected: downloads two files, then prints `<N> simulated genes, <k> up, <k> down` with N ≥ 20 and k ≥ 5. If k < 5, widen the slice (`END = 23_000_000`) and re-run.

- [ ] **Step 3: Verify outputs**

Run: `ls assets/test/reads/*/ && head -3 assets/test/samplesheet.csv && grep -c '>' assets/test/genome.fa && du -sh assets/test --exclude=_downloads`
Expected: six sample folders (ctrl_2 and trt_2 hold only `R1.fastq.gz`); header `sample,fastq_1,fastq_2,condition,batch`; one FASTA record; total under 15 MB.

- [ ] **Step 4: Commit**

```bash
git add assets/test/make_test_data.py assets/test/genome.fa assets/test/genes.gtf assets/test/samplesheet.csv assets/test/de_truth.tsv assets/test/reads
git commit -m "test: add simulated chr22-slice test dataset and its generator

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Config and input validation

**Files:**
- Create: `nextflow.config`, `main.nf`
- Test: `tests/validation.sh`

**Interfaces:**
- Consumes: `assets/test/samplesheet.csv` and reads from Task 3.
- Produces (used by Tasks 5 and 7):
  - Functions in `main.nf`: `checkParams()`, `readSamplesheet(String path) -> List<Map>` (keys trimmed, BOM removed, values trimmed), `resolvePath(String p, base) -> Path` (relative paths resolve against the samplesheet folder), `checkDesign(List rows)`, `checkContrasts(List levels)`, `genomeCache() -> String`, `ensemblUrls(String genome, release) -> Map[fasta, gtf]`.
  - Channel `reads`: `[ [id: String, single_end: Boolean], [Path, Path?] ]`.
  - Params and profiles exactly as in `nextflow.config` below; `params.containers` map keys `python, fastqc, fastp, star, subread, multiqc, report`.

- [ ] **Step 1: Write the failing validation tests**

`tests/validation.sh`:

```bash
#!/usr/bin/env bash
# Input-validation tests. They need no Docker: every check in main.nf runs before any
# process starts, and passing cases use -preview (build the workflow, run nothing).
# Usage: bash tests/validation.sh
set -u
REPO="$(cd "$(dirname "$0")/.." && pwd)"
READS="$REPO/assets/test/reads"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$TMP"
pass=0; fail=0

run() { nextflow run "$REPO/main.nf" -ansi-log false -w "$TMP/work" "$@" 2>&1; }

expect_error() {  # expect_error "<message fragment>" <nextflow args...>
    local msg="$1"; shift
    local out code
    out="$(run "$@")"; code=$?
    if [[ $code -ne 0 && "$out" == *"$msg"* ]]; then echo "PASS  error: $msg"; pass=$((pass + 1))
    else echo "FAIL  expected error '$msg' (exit $code)"; echo "$out" | tail -5; fail=$((fail + 1)); fi
}

expect_ok() {  # expect_ok "<label>" <nextflow args...>
    local label="$1"; shift
    local out code
    out="$(run "$@" -preview)"; code=$?
    if [[ $code -eq 0 ]]; then echo "PASS  ok: $label"; pass=$((pass + 1))
    else echo "FAIL  expected success: $label"; echo "$out" | tail -5; fail=$((fail + 1)); fi
}

sheet() {  # sheet <name> <csv body...>  -> writes $TMP/<name>.csv
    local name="$1"; shift
    printf '%s\n' "$@" > "$TMP/$name.csv"
}

HEADER="sample,fastq_1,fastq_2,condition,batch"
OK1="ctrl_1,$READS/ctrl_1/R1.fastq.gz,$READS/ctrl_1/R2.fastq.gz,control,A"
OK2="trt_1,$READS/trt_1/R1.fastq.gz,$READS/trt_1/R2.fastq.gz,treated,A"
OK3="ctrl_2,$READS/ctrl_2/R1.fastq.gz,,control,B"
OK4="trt_2,$READS/trt_2/R1.fastq.gz,,treated,B"

sheet ok "$HEADER" "$OK1" "$OK2" "$OK3" "$OK4"
sheet no_condition "sample,fastq_1,fastq_2" "ctrl_1,$READS/ctrl_1/R1.fastq.gz,"
sheet duplicate "$HEADER" "$OK1" "$OK1" "$OK2"
sheet missing_fastq "$HEADER" "$OK1" "trt_1,$READS/nope/R1.fastq.gz,,treated,A"
sheet plain_fastq "$HEADER" "$OK1" "trt_1,$READS/trt_1/R1.fastq,,treated,A"
sheet one_condition "$HEADER" "$OK1" "$OK3"
sheet bad_level "$HEADER" "$OK1" "trt_1,$READS/trt_1/R1.fastq.gz,,KO-1,A"
sheet empty_batch "$HEADER" "$OK1" "$OK2" "ctrl_2,$READS/ctrl_2/R1.fastq.gz,,control,"
# Excel-style: UTF-8 BOM + CRLF line endings
printf '\xef\xbb\xbf%s\r\n%s\r\n%s\r\n' "$HEADER" "$OK1" "$OK2" > "$TMP/bom_crlf.csv"
# Relative paths resolve against the samplesheet's own folder
mkdir -p "$TMP/rel" && ln -s "$READS" "$TMP/rel/reads"
sheet rel/relative "$HEADER" "ctrl_1,reads/ctrl_1/R1.fastq.gz,,control,A" "trt_1,reads/trt_1/R1.fastq.gz,,treated,A"

expect_error "Please provide a samplesheet"                 
expect_error "Samplesheet not found"                        --input "$TMP/nope.csv"
expect_error "missing required column(s): condition"        --input "$TMP/no_condition.csv"
expect_error "Duplicate sample name(s) in samplesheet: ctrl_1" --input "$TMP/duplicate.csv"
expect_error "FASTQ file not found"                         --input "$TMP/missing_fastq.csv"
expect_error "must be gzipped"                              --input "$TMP/plain_fastq.csv"
expect_error "at least two conditions"                      --input "$TMP/one_condition.csv"
expect_error "condition values must start with a letter"    --input "$TMP/bad_level.csv"
expect_error "empty design values for sample(s): ctrl_2"    --input "$TMP/empty_batch.csv" --design '~ batch + condition'
expect_error "--strandedness must be one of"                --input "$TMP/ok.csv" --strandedness revrse
expect_error "--genome must be one of"                      --input "$TMP/ok.csv" --genome hg19
expect_error "--fasta and --gtf must be given together"     --input "$TMP/ok.csv" --fasta "$REPO/assets/test/genome.fa"
expect_error "column(s) not in the samplesheet: batch2"     --input "$TMP/ok.csv" --design '~ batch2 + condition'
expect_error "--design must end with 'condition'"           --input "$TMP/ok.csv" --design '~ condition + batch'
expect_error "interaction terms"                            --input "$TMP/ok.csv" --design '~ batch * condition'
expect_error "Contrast 'treated_vs_ctrl' is not valid"      --input "$TMP/ok.csv" --contrasts treated_vs_ctrl
expect_error "--padj_cutoff must be between 0 and 1"        --input "$TMP/ok.csv" --padj_cutoff 5

expect_ok "valid samplesheet"            --input "$TMP/ok.csv" --design '~ batch + condition' --contrasts treated_vs_control
expect_ok "Excel BOM + CRLF samplesheet" --input "$TMP/bom_crlf.csv" --design '~ batch + condition'
expect_ok "relative FASTQ paths"         --input "$TMP/rel/relative.csv"
expect_ok "test profile"                 -profile test
expect_ok "mouse genome"                 --input "$TMP/ok.csv" --genome GRCm39

echo "---- $pass passed, $fail failed"
[[ $fail -eq 0 ]]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `bash tests/validation.sh`
Expected: FAIL lines (no `main.nf` yet), final line `---- 0 passed, 22 failed`.

- [ ] **Step 3: Write `nextflow.config`**

```groovy
// bulk-rnaseq-lite-nf configuration.
// Every value under `params` can be changed on the command line (--name value)
// or in a YAML file passed with -params-file. The command line wins.

manifest {
    name            = 'bulk-rnaseq-lite-nf'
    description     = 'Compact bulk RNA-seq pipeline: FASTQ -> STAR -> featureCounts -> DESeq2/edgeR report'
    mainScript      = 'main.nf'
    nextflowVersion = '>=25.04'
    version         = '0.1.0'
}

params {
    // Input / output
    input           = null            // samplesheet.csv (required)
    outdir          = 'results'

    // Reference genome
    genome          = 'GRCh38'        // GRCh38 (human) or GRCm39 (mouse)
    ensembl_release = 112
    fasta           = null            // custom genome FASTA (needs --gtf too)
    gtf             = null            // custom annotation GTF
    star_index      = null            // existing STAR index folder (built with a GTF)
    genome_cache    = null            // default: ~/.rnaseq-refs
    low_memory      = false           // set by -profile low_memory

    // Reads / alignment
    strandedness    = 'auto'          // auto | forward | reverse | unstranded
    read_length     = 100             // STAR --sjdbOverhang = read_length - 1
    skip_trimming   = false
    save_bam        = false           // copy BAM files to results/star/

    // Differential expression
    design          = '~ condition'   // e.g. '~ batch + condition'
    contrasts       = null            // e.g. 'treated_vs_control,KO_vs_WT'; null = every level vs the first
    padj_cutoff     = 0.05
    lfc_cutoff      = 0.58            // log2 fold change (0.58 ~ 1.5-fold)
    skip_enrichment = false

    // Resource caps: no task asks for more than this
    max_cpus        = 8
    max_memory      = '32.GB'

    // Container images (one pinned image per tool)
    containers = [
        python : 'python:3.12-slim',
        fastqc : 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0',
        fastp  : 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0',
        star   : 'quay.io/biocontainers/star:2.7.11b--h43eeafb_0',
        subread: 'quay.io/biocontainers/subread:2.0.6--he4a0461_0',
        multiqc: 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0',
        report : 'bulk-rnaseq-lite-report:1.0'   // build: docker build -t bulk-rnaseq-lite-report:1.0 docker/
    ]
}

process {
    shell          = ['/bin/bash', '-euo', 'pipefail']
    resourceLimits = [cpus: params.max_cpus, memory: params.max_memory]
    cpus   = 1
    memory = 2.GB

    withLabel: 'process_low'    { cpus = 2; memory = 4.GB }
    withLabel: 'process_medium' { cpus = 4; memory = 8.GB }

    // STAR needs ~32 GB for human/mouse. Out-of-memory kills (exit 137/140) retry once with double memory.
    withName: 'STAR_INDEX|STAR_ALIGN' {
        cpus          = 8
        memory        = { 32.GB * task.attempt }
        errorStrategy = { task.exitStatus in [137, 140] ? 'retry' : 'terminate' }
        maxRetries    = 1
    }
    withName: 'FEATURECOUNTS' {
        cpus          = 4
        memory        = { 4.GB * task.attempt }
        errorStrategy = { task.exitStatus in [137, 140] ? 'retry' : 'terminate' }
        maxRetries    = 1
    }

    withName: 'DOWNLOAD_REFS|INFER_STRAND|MERGE_COUNTS|GENE_NAMES|QC_SUMMARY' { container = params.containers.python }
    withName: 'FASTQC'        { container = params.containers.fastqc }
    withName: 'FASTP'         { container = params.containers.fastp }
    withName: 'STAR_INDEX|STAR_ALIGN' { container = params.containers.star }
    withName: 'FEATURECOUNTS' { container = params.containers.subread }
    withName: 'MULTIQC'       { container = params.containers.multiqc }
    withName: 'REPORT'        { container = params.containers.report }
}

profiles {
    docker {
        docker.enabled    = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }
    // For machines with ~16 GB RAM: sparse STAR index (slower, much less memory).
    low_memory {
        params.low_memory = true
        process {
            withName: 'STAR_INDEX|STAR_ALIGN' { memory = { 14.GB * task.attempt } }
        }
    }
    // Tiny bundled dataset: nextflow run main.nf -profile test,docker
    test {
        params {
            input       = "${projectDir}/assets/test/samplesheet.csv"
            fasta       = "${projectDir}/assets/test/genome.fa"
            gtf         = "${projectDir}/assets/test/genes.gtf"
            genome      = 'GRCh38'
            read_length = 50
            design      = '~ batch + condition'
        }
        process.resourceLimits = [cpus: 2, memory: 6.GB]
    }
}

timeline { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/timeline.html" }
report   { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_report.html" }
trace    { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/trace.txt" }
```

- [ ] **Step 4: Write `main.nf` (checks + reads channel)**

```groovy
#!/usr/bin/env nextflow
/*
 * bulk-rnaseq-lite-nf
 * FASTQ -> FastQC + fastp -> STAR -> featureCounts -> DESeq2 + edgeR + GO -> HTML report
 *
 * Read this file top to bottom: input checks first, then the workflow that wires the
 * modules together. Each module in modules/ is one tool.
 */

// ---------------------------------------------------------------------------
// Input checks. They run before any process starts, so mistakes fail in seconds.
// ---------------------------------------------------------------------------

def validStrandedness() {
    return ['auto', 'forward', 'reverse', 'unstranded']
}

def supportedGenomes() {
    return ['GRCh38', 'GRCm39']
}

def checkParams() {
    if (!params.input) {
        error("Please provide a samplesheet with --input samplesheet.csv")
    }
    if (!(params.strandedness in validStrandedness())) {
        error("--strandedness must be one of ${validStrandedness().join(', ')} (got '${params.strandedness}')")
    }
    if (params.fasta && !params.gtf) {
        error("--fasta and --gtf must be given together")
    }
    if (params.gtf && !params.fasta && !params.star_index) {
        error("--gtf needs --fasta (to build an index) or --star_index")
    }
    if (!params.gtf && !(params.genome in supportedGenomes())) {
        error("--genome must be one of ${supportedGenomes().join(', ')} (got '${params.genome}'). For other species give --fasta and --gtf.")
    }
    if (!(params.padj_cutoff > 0 && params.padj_cutoff <= 1)) {
        error("--padj_cutoff must be between 0 and 1 (got ${params.padj_cutoff})")
    }
    if (params.lfc_cutoff < 0) {
        error("--lfc_cutoff must be 0 or more (got ${params.lfc_cutoff})")
    }
}

// Relative FASTQ paths are read relative to the samplesheet's folder.
def resolvePath(String p, base) {
    return (p.startsWith('/') || p.contains('://')) ? file(p) : base.resolve(p)
}

def readSamplesheet(String path) {
    def sheet = file(path)
    if (!sheet.exists()) {
        error("Samplesheet not found: ${path}")
    }
    // Trim keys and values; drop the invisible BOM that Excel adds to the first header.
    def rows = sheet.splitCsv(header: true).collect { row ->
        row.collectEntries { k, v -> [k.replace('﻿', '').trim(), (v ?: '').trim()] }
    }
    if (!rows) {
        error("Samplesheet ${path} has no samples")
    }
    def missing = ['sample', 'fastq_1', 'condition'].findAll { c -> !rows[0].containsKey(c) }
    if (missing) {
        error("Samplesheet is missing required column(s): ${missing.join(', ')}")
    }
    def names = rows.collect { r -> r.sample }
    def dups = names.findAll { n -> names.count(n) > 1 }.unique()
    if (dups) {
        error("Duplicate sample name(s) in samplesheet: ${dups.join(', ')}")
    }
    rows.each { r ->
        if (!r.sample.matches('[A-Za-z0-9_.\\-]+')) {
            error("Sample name '${r.sample}' may only contain letters, digits, '_', '.' and '-'")
        }
        if (!r.condition.matches('[A-Za-z][A-Za-z0-9_.]*') || r.condition.contains('_vs_')) {
            error("condition values must start with a letter, contain only letters, digits, '_' or '.', and not contain '_vs_' (got '${r.condition}' for sample ${r.sample})")
        }
        [r.fastq_1, r.fastq_2].findAll { p -> p }.each { p ->
            if (!(p.endsWith('.fastq.gz') || p.endsWith('.fq.gz'))) {
                error("FASTQ files must be gzipped and end in .fastq.gz or .fq.gz: ${p}")
            }
            if (!resolvePath(p, sheet.parent).exists()) {
                error("FASTQ file not found: ${p} (sample ${r.sample})")
            }
        }
    }
    return rows
}

def checkDesign(List rows) {
    def design = params.design.toString().trim()
    if (!design.startsWith('~')) {
        error("--design must start with '~', e.g. '~ batch + condition' (got '${design}')")
    }
    if (design.contains('*') || design.contains(':')) {
        error("--design interaction terms (* or :) are not supported (got '${design}')")
    }
    def vars = design.substring(1).tokenize('+').collect { v -> v.trim() }.findAll { v -> v }
    if (!vars || vars.last() != 'condition') {
        error("--design must end with 'condition', e.g. '~ batch + condition' (got '${design}')")
    }
    def missing = vars.findAll { v -> !rows[0].containsKey(v) }
    if (missing) {
        error("--design uses column(s) not in the samplesheet: ${missing.join(', ')}")
    }
    def empty = rows.findAll { r -> vars.any { v -> !r[v] } }.collect { r -> r.sample }
    if (empty) {
        error("Samplesheet has empty design values for sample(s): ${empty.join(', ')}")
    }
}

def checkContrasts(List levels) {
    if (levels.size() < 2) {
        error("The samplesheet needs at least two conditions (found: ${levels.join(', ')})")
    }
    if (!params.contrasts) {
        return
    }
    params.contrasts.toString().tokenize(',').collect { c -> c.trim() }.each { c ->
        def parts = c.split('_vs_')
        if (parts.size() != 2 || !(parts[0] in levels) || !(parts[1] in levels)) {
            error("Contrast '${c}' is not valid. Write LEVEL_vs_LEVEL using: ${levels.join(', ')}")
        }
    }
}

// ---------------------------------------------------------------------------
// Reference helpers
// ---------------------------------------------------------------------------

def genomeCache() {
    return params.genome_cache ?: "${System.getenv('HOME')}/.rnaseq-refs"
}

def ensemblUrls(String genome, release) {
    def species = [GRCh38: 'homo_sapiens', GRCm39: 'mus_musculus'][genome]
    def prefix  = [GRCh38: 'Homo_sapiens', GRCm39: 'Mus_musculus'][genome]
    def base    = "https://ftp.ensembl.org/pub/release-${release}"
    return [
        fasta: "${base}/fasta/${species}/dna/${prefix}.${genome}.dna.primary_assembly.fa.gz",
        gtf  : "${base}/gtf/${species}/${prefix}.${genome}.${release}.gtf.gz"
    ]
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
    checkParams()
    def rows = readSamplesheet(params.input)
    checkDesign(rows)
    checkContrasts(rows.collect { r -> r.condition }.unique())

    // One item per sample: [ [id, single_end], [fastq_1, (fastq_2)] ]
    def sheet_dir = file(params.input).parent
    def reads = channel.fromList(rows).map { r ->
        def files = [r.fastq_1, r.fastq_2].findAll { p -> p }.collect { p -> resolvePath(p, sheet_dir) }
        [[id: r.sample, single_end: files.size() == 1], files]
    }
}
```

- [ ] **Step 5: Lint and run the validation tests**

Run: `nextflow lint main.nf nextflow.config && bash tests/validation.sh`
Expected: lint reports no errors; `---- 22 passed, 0 failed`.
If lint rejects a construct, rewrite it in the form lint suggests and re-run (don't disable lint).

- [ ] **Step 6: Commit**

```bash
git add nextflow.config main.nf tests/validation.sh
git commit -m "feat: add config, profiles and fail-fast input validation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Upstream modules (references → MultiQC) and wiring

**Files:**
- Create: `modules/download_refs.nf`, `modules/star_index.nf`, `modules/fastqc.nf`, `modules/fastp.nf`, `modules/star_align.nf`, `modules/infer_strand.nf`, `modules/featurecounts.nf`, `modules/helpers.nf`, `modules/multiqc.nf`, `assets/NO_FILE`
- Modify: `main.nf` (add includes at top; extend the `workflow` block)

**Interfaces:**
- Consumes: `reads` channel, `genomeCache()`, `ensemblUrls()` from Task 4; `bin/*.py` from Task 2.
- Produces (published, used by Task 7 and the Rmd defaults):
  - `results/featurecounts/counts.tsv`, `results/featurecounts/gene_names.tsv`
  - `results/qc/qc_summary.tsv`
  - `results/star/strandedness.txt`, `results/star/strandedness_mqc.tsv` (auto mode only)
  - `results/multiqc/multiqc_report.html`
  - `results/pipeline_info/samplesheet.csv`, `results/pipeline_info/versions.tsv`
  - Workflow variables for Task 7: `MERGE_COUNTS.out`, `GENE_NAMES.out`, `QC_SUMMARY.out`, `strandedness` (value channel of String), `strand_report` (value channel of Path; `assets/NO_FILE` when strandedness was set by the user), `samplesheet_copy`, `versions`.

- [ ] **Step 1: Write `modules/download_refs.nf`**

```groovy
// Download the Ensembl genome FASTA and GTF once. storeDir makes Nextflow skip this
// process entirely when both files already exist in the cache folder.
process DOWNLOAD_REFS {
    tag "${store_dir}"
    storeDir "${store_dir}"

    input:
    val fasta_url
    val gtf_url
    val store_dir

    output:
    path 'genome.fa', emit: fasta
    path 'genes.gtf', emit: gtf

    script:
    """
    python3 -c "import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], 'genome.fa.gz')" ${fasta_url}
    python3 -c "import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], 'genes.gtf.gz')" ${gtf_url}
    gunzip genome.fa.gz genes.gtf.gz
    """
}
```

- [ ] **Step 2: Write `modules/star_index.nf`**

```groovy
// Build the STAR genome index once and keep it in the cache (storeDir).
// -profile low_memory builds a sparse index that needs ~16 GB instead of ~32 GB.
process STAR_INDEX {
    tag "${params.low_memory ? 'sparse' : 'full'} index"
    storeDir "${store_dir}"

    input:
    path fasta
    path gtf
    val store_dir

    output:
    path "${params.low_memory ? 'star_index_sparse3' : 'star_index'}", emit: index

    script:
    def index_dir = params.low_memory ? 'star_index_sparse3' : 'star_index'
    def sparse = params.low_memory ? '--genomeSAsparseD 3' : ''
    """
    # STAR cannot read gzipped references
    if [[ ${fasta} == *.gz ]]; then gunzip -c ${fasta} > ref.fa; else ln -s ${fasta} ref.fa; fi
    if [[ ${gtf} == *.gz ]]; then gunzip -c ${gtf} > ref.gtf; else ln -s ${gtf} ref.gtf; fi

    # Small genomes need a smaller suffix-array index: min(14, log2(genome length)/2 - 1)
    genome_len=\$(grep -v '>' ref.fa | tr -d '\\n' | wc -c)
    sa_bases=\$(awk -v len="\$genome_len" 'BEGIN { n = int(log(len) / log(2) / 2 - 1); print (n < 14 ? n : 14) }')

    mkdir ${index_dir}
    STAR --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index_dir} \\
        --genomeFastaFiles ref.fa \\
        --sjdbGTFfile ref.gtf \\
        --sjdbOverhang ${params.read_length - 1} \\
        --genomeSAindexNbases \$sa_bases \\
        --limitGenomeGenerateRAM ${task.memory.toBytes()} \\
        ${sparse}
    """
}
```

- [ ] **Step 3: Write `modules/fastqc.nf`**

```groovy
// Quality report for the raw reads.
// Inputs are staged into input1/, input2/ so two samples with the same FASTQ file name
// never collide; the links give each report a sample-based name.
process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads, stageAs: 'input?/*')

    output:
    path '*_fastqc.zip', emit: zip
    path '*_fastqc.html', emit: html

    script:
    def r = [reads].flatten()   // always a list, for single- and paired-end
    def links = meta.single_end
        ? "ln -s ${r[0]} ${meta.id}.fastq.gz"
        : "ln -s ${r[0]} ${meta.id}_1.fastq.gz && ln -s ${r[1]} ${meta.id}_2.fastq.gz"
    """
    ${links}
    fastqc --threads ${task.cpus} --quiet ${meta.id}*.fastq.gz
    """
}
```

- [ ] **Step 4: Write `modules/fastp.nf`**

```groovy
// Adapter and quality trimming. fastp detects adapters by itself and writes a JSON
// report that MultiQC and the QC summary read.
process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(meta), path(reads, stageAs: 'input?/*')

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html", emit: html

    script:
    def r = [reads].flatten()
    def io = meta.single_end
        ? "--in1 ${r[0]} --out1 ${meta.id}.trim.fastq.gz"
        : "--in1 ${r[0]} --in2 ${r[1]} --out1 ${meta.id}_1.trim.fastq.gz --out2 ${meta.id}_2.trim.fastq.gz --detect_adapter_for_pe"
    """
    fastp ${io} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html
    """
}
```

- [ ] **Step 5: Write `modules/star_align.nf`**

```groovy
// Align reads with STAR. --quantMode GeneCounts also writes ReadsPerGene.out.tab,
// which INFER_STRAND uses to detect the library strandedness.
process STAR_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/star", mode: 'copy',
        saveAs: { fn -> (fn.endsWith('.bam') && !params.save_bam) ? null : fn }

    input:
    tuple val(meta), path(reads, stageAs: 'input?/*')
    path index

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${meta.id}.Log.final.out", emit: log
    path "${meta.id}.ReadsPerGene.out.tab", emit: gene_counts

    script:
    def files = [reads].flatten().join(' ')
    """
    STAR --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${files} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS nM \\
        --quantMode GeneCounts \\
        --limitBAMsortRAM ${task.memory.toBytes().intdiv(2)}
    """
}
```

- [ ] **Step 6: Write `modules/infer_strand.nf`**

```groovy
// Pick one strandedness for the whole run from all samples' STAR gene counts.
// Stops the run if samples disagree or the call is ambiguous (see bin/infer_strand.py).
process INFER_STRAND {
    label 'process_low'
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    path gene_counts

    output:
    path 'strandedness.txt', emit: strand
    path 'strandedness_mqc.tsv', emit: mqc

    script:
    """
    infer_strand.py ${gene_counts}
    """
}
```

- [ ] **Step 7: Write `modules/featurecounts.nf`**

```groovy
// Count reads per gene (exons summarised by gene_id). Runs per sample because one
// featureCounts call cannot mix paired-end (-p) and single-end BAMs.
process FEATURECOUNTS {
    tag "${meta.id}"
    publishDir "${params.outdir}/featurecounts/per_sample", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path gtf
    val strandedness

    output:
    path "${meta.id}.featureCounts.txt", emit: counts
    path "${meta.id}.featureCounts.txt.summary", emit: summary

    script:
    def strand_flag = [unstranded: 0, forward: 1, reverse: 2][strandedness]
    def paired = meta.single_end ? '' : '-p --countReadPairs'
    """
    featureCounts -T ${task.cpus} \\
        -a ${gtf} -t exon -g gene_id \\
        -s ${strand_flag} ${paired} \\
        -o ${meta.id}.featureCounts.txt \\
        ${bam}
    """
}
```

- [ ] **Step 8: Write `modules/helpers.nf`**

```groovy
// Small steps that run the Python helpers in bin/.

// All per-sample featureCounts files -> one gene x sample matrix (counts.tsv).
process MERGE_COUNTS {
    label 'process_low'
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    path count_files

    output:
    path 'counts.tsv'

    script:
    """
    merge_counts.py ${count_files}
    """
}

// gene_id -> gene symbol and biotype, read straight from the GTF.
process GENE_NAMES {
    label 'process_low'
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    path gtf

    output:
    path 'gene_names.tsv'

    script:
    """
    gene_names.py ${gtf}
    """
}

// One QC table (reads, % trimmed, % mapped, % assigned) for the report.
process QC_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path qc_files

    output:
    path 'qc_summary.tsv'

    script:
    """
    qc_summary.py ${qc_files}
    """
}
```

- [ ] **Step 9: Write `modules/multiqc.nf`**

```groovy
// One HTML page with FastQC, fastp, STAR, featureCounts and strandedness results.
process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc . --force --filename multiqc_report.html
    """
}
```

- [ ] **Step 10: Create the placeholder file**

Run: `touch assets/NO_FILE`

- [ ] **Step 11: Add includes to `main.nf`**

Insert directly after the header comment block (before `// Input checks`):

```groovy
include { DOWNLOAD_REFS } from './modules/download_refs'
include { STAR_INDEX    } from './modules/star_index'
include { FASTQC        } from './modules/fastqc'
include { FASTP         } from './modules/fastp'
include { STAR_ALIGN    } from './modules/star_align'
include { INFER_STRAND  } from './modules/infer_strand'
include { FEATURECOUNTS } from './modules/featurecounts'
include { MERGE_COUNTS; GENE_NAMES; QC_SUMMARY } from './modules/helpers'
include { MULTIQC       } from './modules/multiqc'
```

- [ ] **Step 12: Extend the `workflow` block in `main.nf`**

Append inside `workflow { ... }`, after the `reads = ...` statement:

```groovy
    // ---- Reference genome ------------------------------------------------
    // Ensembl downloads live in the shared cache; indexes for custom FASTA/GTF live in <outdir>/reference.
    def cache_dir = params.gtf
        ? "${params.outdir}/reference"
        : "${genomeCache()}/${params.genome}/ensembl_${params.ensembl_release}"
    def fasta = channel.empty()
    def gtf = channel.empty()
    if (params.gtf) {
        gtf = channel.value(file(params.gtf, checkIfExists: true))
        if (params.fasta) {
            fasta = channel.value(file(params.fasta, checkIfExists: true))
        }
    } else {
        def urls = ensemblUrls(params.genome, params.ensembl_release)
        DOWNLOAD_REFS(urls.fasta, urls.gtf, cache_dir)
        fasta = DOWNLOAD_REFS.out.fasta
        gtf = DOWNLOAD_REFS.out.gtf
    }
    def index = channel.empty()
    if (params.star_index) {
        index = channel.value(file(params.star_index, checkIfExists: true))
    } else {
        STAR_INDEX(fasta, gtf, cache_dir)
        index = STAR_INDEX.out.index
    }

    // ---- QC, trimming, alignment ----------------------------------------
    FASTQC(reads)
    def trimmed = reads
    def fastp_json = channel.empty()
    if (!params.skip_trimming) {
        FASTP(reads)
        trimmed = FASTP.out.reads
        fastp_json = FASTP.out.json
    }
    STAR_ALIGN(trimmed, index)

    // ---- Strandedness: auto-detect from STAR gene counts, or use --strandedness
    def strandedness = channel.value(params.strandedness)
    def strand_report = channel.value(file("${projectDir}/assets/NO_FILE"))
    if (params.strandedness == 'auto') {
        INFER_STRAND(STAR_ALIGN.out.gene_counts.collect())
        strandedness = INFER_STRAND.out.strand.map { f -> f.text.trim() }
        strand_report = INFER_STRAND.out.mqc
    }

    // ---- Counting and QC summaries --------------------------------------
    FEATURECOUNTS(STAR_ALIGN.out.bam, gtf, strandedness)
    MERGE_COUNTS(FEATURECOUNTS.out.counts.collect())
    GENE_NAMES(gtf)
    QC_SUMMARY(fastp_json.mix(STAR_ALIGN.out.log, FEATURECOUNTS.out.summary).collect())
    MULTIQC(
        FASTQC.out.zip
            .mix(fastp_json, STAR_ALIGN.out.log, FEATURECOUNTS.out.summary, strand_report)
            .collect()
    )

    // ---- Run records: samplesheet copy + tool versions (pipeline_info/) ---
    def samplesheet_copy = channel.fromPath(params.input)
        .collectFile(name: 'samplesheet.csv', storeDir: "${params.outdir}/pipeline_info")
    def version_lines = ["nextflow\t${workflow.nextflow.version}", "pipeline\t${workflow.manifest.version}"] +
        params.containers.collect { k, v -> "${k}\t${v}" }
    def versions = channel.fromList(version_lines)
        .collectFile(name: 'versions.tsv', newLine: true, sort: false, storeDir: "${params.outdir}/pipeline_info")
```

- [ ] **Step 13: Lint and re-run validation**

Run: `nextflow lint main.nf modules/ && bash tests/validation.sh`
Expected: no lint errors; `---- 22 passed, 0 failed`.

- [ ] **Step 14: Run the test profile end to end**

Run: `nextflow run main.nf -profile test,docker`
Expected: completes; `DOWNLOAD_REFS` does not appear (custom FASTA/GTF); all other processes succeed.

- [ ] **Step 15: Check the outputs**

Run:
```bash
cat results/star/strandedness.txt
head -2 results/featurecounts/counts.tsv
cat results/qc/qc_summary.tsv
ls results/multiqc/multiqc_report.html results/pipeline_info/samplesheet.csv results/pipeline_info/versions.tsv
ls results/fastqc/ | head
```
Expected: `reverse`; header `gene_id	ctrl_1	ctrl_2	ctrl_3	trt_1	trt_2	trt_3`; `pct_uniquely_mapped` > 80 and `pct_assigned` > 50 for every sample; the three files exist; FastQC files named `ctrl_1_1_fastqc.html`, `ctrl_2_fastqc.html` etc. (sample-based, no collisions from the identical `R1.fastq.gz` names).

- [ ] **Step 16: Review Focus 5 — the index is reused on a second run**

Run:
```bash
before=$(stat -c %Y results/reference/star_index/SA)
nextflow run main.nf -profile test,docker
after=$(stat -c %Y results/reference/star_index/SA)
[[ "$before" == "$after" ]] && echo "index reused" || echo "INDEX REBUILT"
```
Expected: `index reused`.

- [ ] **Step 17: Explicit strandedness path**

Run: `nextflow run main.nf -profile test,docker --strandedness reverse --outdir results_explicit && ls results_explicit/star/ | grep -c strandedness; rm -rf results_explicit`
Expected: completes; prints `0` (no strandedness files, since detection was skipped).

- [ ] **Step 18: Commit**

```bash
git add modules/ assets/NO_FILE main.nf
git commit -m "feat: add upstream modules (refs, FastQC, fastp, STAR, strand, featureCounts, MultiQC)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Report Docker image

**Files:**
- Create: `docker/Dockerfile`

**Interfaces:**
- Produces: local image `bulk-rnaseq-lite-report:1.0` (name matches `params.containers.report`) containing R, pandoc and the packages listed below.

- [ ] **Step 1: Write `docker/Dockerfile`**

```dockerfile
# R environment for the bulk-rnaseq-lite-nf HTML report.
# Build once:  docker build -t bulk-rnaseq-lite-report:1.0 docker/
# Bioconductor 3.22 pins package versions; the base image already has R, pandoc and rmarkdown.
FROM bioconductor/bioconductor_docker:RELEASE_3_22

RUN R -q -e 'pkgs <- c("DESeq2", "edgeR", "limma", "apeglm", "clusterProfiler", "enrichplot", \
                       "org.Hs.eg.db", "org.Mm.eg.db", "rmarkdown", "knitr", "ggplot2", "pheatmap", "DT"); \
             BiocManager::install(pkgs, ask = FALSE, update = FALSE); \
             missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
             if (length(missing) > 0) stop("Failed to install: ", paste(missing, collapse = ", "))'
```

- [ ] **Step 2: Build the image**

Run: `docker build -t bulk-rnaseq-lite-report:1.0 docker/`
Expected: finishes without `Failed to install`. (10–30 min the first time.)

- [ ] **Step 3: Verify the packages load and pandoc exists**

Run: `docker run --rm bulk-rnaseq-lite-report:1.0 Rscript -e 'suppressPackageStartupMessages({library(DESeq2); library(edgeR); library(apeglm); library(clusterProfiler); library(org.Hs.eg.db); library(org.Mm.eg.db)}); cat(rmarkdown::pandoc_available(), "\n")'`
Expected: `TRUE`.

- [ ] **Step 4: Commit**

```bash
git add docker/Dockerfile
git commit -m "feat: add Docker image for the R report

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: R Markdown report, REPORT module and end-to-end check

**Files:**
- Create: `report/rnaseq_report.Rmd`, `report/_contrast.Rmd`, `modules/report.nf`, `tests/check_results.py`
- Modify: `main.nf` (include REPORT; call it at the end of the workflow)

**Interfaces:**
- Consumes: Task 5 outputs and workflow variables (`MERGE_COUNTS.out`, `GENE_NAMES.out`, `QC_SUMMARY.out`, `strand_report`, `strandedness`, `samplesheet_copy`, `versions`); image from Task 6.
- Produces: `results/report/report.html`, `results/report/tables/<contrast>_{deseq2,edger}.csv`, and `<contrast>_go_{up,down}.csv` when GO terms are found. DESeq2 CSV columns: `gene_id, gene_name, baseMean, log2FC, log2FC_shrunk, pvalue, padj, significant`. edgeR CSV columns: `gene_id, gene_name, logCPM, log2FC, pvalue, padj, significant`.
- Rmd params (names fixed; pipeline and standalone use the same ones): `counts, samplesheet, gene_names, qc_summary, strand_report, versions, multiqc_link, design, contrasts, padj_cutoff, lfc_cutoff, genome, ensembl_release, strandedness, skip_enrichment, tables_dir`.

- [ ] **Step 1: Write the failing end-to-end check `tests/check_results.py`**

```python
#!/usr/bin/env python3
"""End-to-end check after `nextflow run main.nf -profile test,docker`.

The test data were simulated with known up/down genes (assets/test/de_truth.tsv);
both DESeq2 and edgeR should recover most of them with few false positives.
Usage: python3 tests/check_results.py [results_dir]
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
results = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results"
MAX_FALSE_POSITIVES = 3

truth = {}
with open(ROOT / "assets/test/de_truth.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        truth[row["gene_id"]] = row["direction"]

strand = (results / "star/strandedness.txt").read_text().strip()
assert strand == "reverse", f"expected reverse strandedness, got {strand}"

for method in ("deseq2", "edger"):
    with open(results / f"report/tables/treated_vs_control_{method}.csv") as fh:
        rows = list(csv.DictReader(fh))
    found = {r["gene_id"]: ("up" if float(r["log2FC"]) > 0 else "down")
             for r in rows if r["significant"] == "TRUE"}
    for direction in ("up", "down"):
        expected = [g for g, d in truth.items() if d == direction]
        hits = sum(found.get(g) == direction for g in expected)
        assert hits >= len(expected) / 2, f"{method}: only {hits}/{len(expected)} simulated {direction} genes found"
    false_pos = [g for g in found if g not in truth]
    assert len(false_pos) <= MAX_FALSE_POSITIVES, f"{method}: too many false positives: {false_pos}"
    print(f"ok {method}: {len(found)} significant, {len(false_pos)} false positives")

size = (results / "report/report.html").stat().st_size
assert size > 100_000, f"report.html is suspiciously small ({size} bytes)"
print(f"ok report.html ({size // 1024} KB), strandedness={strand}")
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python3 tests/check_results.py`
Expected: FAIL with `FileNotFoundError` for `report/tables/treated_vs_control_deseq2.csv`.

- [ ] **Step 3: Write `report/rnaseq_report.Rmd`**

````markdown
---
title: "Bulk RNA-seq report"
date: "`r format(Sys.time(), '%Y-%m-%d %H:%M')`"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 2
    code_folding: hide
    self_contained: true
params:
  counts: "../results/featurecounts/counts.tsv"
  samplesheet: "../results/pipeline_info/samplesheet.csv"
  gene_names: "../results/featurecounts/gene_names.tsv"
  qc_summary: "../results/qc/qc_summary.tsv"
  strand_report: "../results/star/strandedness_mqc.tsv"
  versions: "../results/pipeline_info/versions.tsv"
  multiqc_link: "../results/multiqc/multiqc_report.html"
  design: "~ condition"
  contrasts: ""
  padj_cutoff: 0.05
  lfc_cutoff: 0.58
  genome: "GRCh38"
  ensembl_release: "112"
  strandedness: "auto"
  skip_enrichment: false
  tables_dir: "tables"
---

```{r setup, include = FALSE}
# HOW TO RUN THIS REPORT WITHOUT NEXTFLOW
#   1. Open this file in RStudio (or use the Docker image, see docs/usage.md).
#   2. Check the paths under `params:` above point at your pipeline results folder.
#      The defaults work when results/ sits next to the report/ folder.
#   3. Change what you need (contrasts, cutoffs, design) and click "Knit".
#   From R instead:
#     rmarkdown::render("rnaseq_report.Rmd", params = list(contrasts = "treated_vs_control"))
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 4.5)
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(DT)
})
theme_set(theme_bw())

# Thresholds used for QC flags and GO. Change them here if your data need it.
MIN_PCT_UNIQUE   <- 60  # % uniquely mapped reads below this is flagged
MIN_PCT_ASSIGNED <- 50  # % reads assigned to genes below this is flagged
MIN_GO_GENES     <- 5   # fewer significant genes than this: skip GO for that list
REPORT_DIR <- dirname(knitr::current_input(dir = TRUE))
```

This report was made by **bulk-rnaseq-lite-nf**. Click a **Code** button to see the R code behind any result.

# Warnings and checks

```{r load-data}
# Read the count matrix: one row per gene, one column per sample.
counts_df <- read.delim(params$counts, check.names = FALSE)
counts <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts) <- counts_df$gene_id

# Read the samplesheet and put its rows in the same order as the count columns.
samples <- read.csv(params$samplesheet, stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8-BOM", strip.white = TRUE)
missing_samples <- setdiff(colnames(counts), samples$sample)
if (length(missing_samples) > 0) {
  stop("These count columns are not in the samplesheet: ", paste(missing_samples, collapse = ", "))
}
# Condition levels keep the samplesheet order; the first one is the reference.
condition_levels <- unique(samples$condition)
samples <- samples[match(colnames(counts), samples$sample), , drop = FALSE]
rownames(samples) <- samples$sample
samples$condition <- factor(samples$condition,
                            levels = condition_levels[condition_levels %in% samples$condition])

# Gene symbols (from the GTF), QC numbers and tool versions written by the pipeline.
gene_names <- read.delim(params$gene_names)
symbol_of <- setNames(gene_names$gene_name, gene_names$gene_id)
qc <- read.delim(params$qc_summary, check.names = FALSE, na.strings = "NA")
versions <- if (file.exists(params$versions)) {
  read.delim(params$versions, header = FALSE, col.names = c("tool", "version"))
}
strand_tab <- if (file.exists(params$strand_report) && basename(params$strand_report) != "NO_FILE") {
  read.delim(params$strand_report, comment.char = "#")
}
strand_used <- if (!is.null(strand_tab)) {
  paste(unique(strand_tab$call), "(auto-detected)")
} else {
  paste(params$strandedness, "(set by user)")
}
```

```{r checks}
report_warnings <- character()

# Condition names must be valid R names so DESeq2/edgeR can name coefficients.
lv <- levels(samples$condition)
if (any(make.names(lv) != lv)) {
  stop("Condition values must start with a letter and use only letters, digits, '.' or '_': ",
       paste(lv, collapse = ", "))
}

# Turn the design text (e.g. "~ batch + condition") into a list of columns.
design_vars <- trimws(strsplit(sub("^\\s*~", "", params$design), "+", fixed = TRUE)[[1]])
design_vars <- design_vars[nzchar(design_vars)]
if (length(design_vars) == 0 || tail(design_vars, 1) != "condition") {
  stop("The design must end with 'condition', e.g. '~ batch + condition'.")
}
not_found <- setdiff(design_vars, colnames(samples))
if (length(not_found) > 0) {
  stop("Design columns missing from the samplesheet: ", paste(not_found, collapse = ", "))
}
# Text columns (like batch) become factors; numeric columns (like RIN) stay numeric covariates.
for (v in setdiff(design_vars, "condition")) {
  if (is.character(samples[[v]])) samples[[v]] <- factor(samples[[v]])
}

# Check 1: can the design be estimated? A batch that lines up exactly with condition
# cannot be separated from it, so we fall back to '~ condition'.
mm <- model.matrix(as.formula(paste("~", paste(design_vars, collapse = " + "))), samples)
if (qr(mm)$rank < ncol(mm)) {
  report_warnings <- c(report_warnings, sprintf(
    "The design '%s' cannot be estimated (a covariate is confounded with condition). It was replaced by '~ condition'.",
    params$design))
  design_vars <- "condition"
  mm <- model.matrix(~ condition, samples)
}
design_formula <- as.formula(paste("~", paste(design_vars, collapse = " + ")))

# Check 2: replicates per group.
reps <- table(samples$condition)
if (any(reps < 2)) {
  report_warnings <- c(report_warnings, sprintf(
    "Groups with fewer than 2 replicates: %s. Their results are unreliable.",
    paste(names(reps)[reps < 2], collapse = ", ")))
}

# Check 3: testing needs more samples than model coefficients.
can_test <- nrow(mm) > ncol(mm)
if (!can_test) {
  report_warnings <- c(report_warnings,
    "Not enough samples to estimate variability (no residual degrees of freedom). Differential expression was skipped.")
}

# Check 4: low assignment usually means wrong strandedness or annotation.
low_assigned <- qc$sample[!is.na(qc$pct_assigned) & qc$pct_assigned < MIN_PCT_ASSIGNED]
if (length(low_assigned) > 0) {
  report_warnings <- c(report_warnings, sprintf(
    "Less than %d%% of reads were assigned to genes in: %s. Check strandedness and annotation.",
    MIN_PCT_ASSIGNED, paste(low_assigned, collapse = ", ")))
}

# Remove genes with too few counts to test. filterByExpr keeps genes with enough
# counts in at least as many samples as the smallest group.
y_all <- DGEList(counts)
keep <- filterByExpr(y_all, group = samples$condition)
counts_f <- counts[keep, , drop = FALSE]
if (sum(keep) < 2) {
  can_test <- FALSE
  report_warnings <- c(report_warnings, "Fewer than 2 genes passed filtering. Differential expression was skipped.")
}

# Contrasts: "A_vs_B" compares A with B (positive log2FC = higher in A).
# Empty = every condition against the first (reference) condition.
contrast_list <- if (is.null(params$contrasts) || !nzchar(params$contrasts)) {
  lapply(lv[-1], function(l) c(l, lv[1]))
} else {
  lapply(trimws(strsplit(params$contrasts, ",")[[1]]), function(x) {
    ab <- strsplit(x, "_vs_", fixed = TRUE)[[1]]
    if (length(ab) != 2 || !all(ab %in% lv)) {
      stop(sprintf("Contrast '%s' is not valid. Write LEVEL_vs_LEVEL using: %s", x, paste(lv, collapse = ", ")))
    }
    ab
  })
}
names(contrast_list) <- vapply(contrast_list, paste, character(1), collapse = "_vs_")
```

```{r show-warnings, echo = FALSE, results = "asis"}
if (length(report_warnings) > 0) {
  cat(paste0("> **Warning:** ", report_warnings, collapse = "\n>\n"), "\n")
} else {
  cat("All checks passed.\n")
}
```

# 1. Run summary

```{r run-summary}
# The settings this report was made with.
run_info <- data.frame(
  setting = c("Genome", "Ensembl release", "Strandedness", "Design (used)", "Contrasts",
              "padj cutoff", "|log2FC| cutoff", "GO enrichment"),
  value = c(params$genome, as.character(params$ensembl_release), strand_used,
            paste(deparse(design_formula), collapse = ""), paste(names(contrast_list), collapse = ", "),
            params$padj_cutoff, params$lfc_cutoff,
            if (isTRUE(params$skip_enrichment)) "skipped" else "GO BP/MF/CC (DESeq2 significant genes)"))
knitr::kable(run_info)
```

```{r versions}
# Tool versions (container images) used by the pipeline.
if (!is.null(versions)) knitr::kable(versions)
```

```{r sample-table}
# The samples and their design columns.
knitr::kable(samples[, setdiff(colnames(samples), c("fastq_1", "fastq_2")), drop = FALSE], row.names = FALSE)
```

# 2. Upstream QC

Detailed per-tool plots are in the [MultiQC report](`r params$multiqc_link`).

```{r qc-table}
# One row per sample. Flagged: uniquely mapped < MIN_PCT_UNIQUE % or assigned < MIN_PCT_ASSIGNED %.
qc$flag <- ifelse((!is.na(qc$pct_uniquely_mapped) & qc$pct_uniquely_mapped < MIN_PCT_UNIQUE) |
                  (!is.na(qc$pct_assigned) & qc$pct_assigned < MIN_PCT_ASSIGNED), "CHECK", "")
knitr::kable(qc, digits = 1, row.names = FALSE)
```

```{r qc-plot, fig.height = 0.4 * nrow(qc) + 2}
# Mapping and gene-assignment rates. Dashed lines are the flag thresholds.
qc_long <- rbind(
  data.frame(sample = qc$sample, metric = "% uniquely mapped", value = qc$pct_uniquely_mapped),
  data.frame(sample = qc$sample, metric = "% assigned to genes", value = qc$pct_assigned))
ggplot(qc_long, aes(sample, value, fill = metric)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = c(MIN_PCT_UNIQUE, MIN_PCT_ASSIGNED), linetype = "dashed") +
  coord_flip() + labs(x = NULL, y = "%", fill = NULL)
```

```{r strand-table}
# Strandedness per sample (auto mode only): forward fraction ~0.5 unstranded, high forward, low reverse.
if (!is.null(strand_tab)) knitr::kable(strand_tab, row.names = FALSE)
```

# 3. Count QC

```{r library-sizes, fig.height = 0.4 * ncol(counts) + 2}
# Library size (reads counted on genes) and number of genes with at least one read.
lib <- data.frame(sample = colnames(counts), million_reads = colSums(counts) / 1e6,
                  genes_detected = colSums(counts > 0), condition = samples$condition)
ggplot(lib, aes(sample, million_reads, fill = condition)) + geom_col() + coord_flip() +
  labs(x = NULL, y = "Million reads assigned to genes")
ggplot(lib, aes(sample, genes_detected, fill = condition)) + geom_col() + coord_flip() +
  labs(x = NULL, y = "Genes with at least one read")
```

```{r logcpm-distribution}
# Distribution of log2 counts-per-million per sample. Boxes should look similar.
logcpm <- cpm(y_all, log = TRUE)
ggplot(data.frame(sample = rep(colnames(logcpm), each = nrow(logcpm)), value = as.vector(logcpm)),
       aes(sample, value)) +
  geom_boxplot(outlier.size = 0.3) + coord_flip() + labs(x = NULL, y = "log2 CPM")
```

**Gene filtering:** `r sum(keep)` of `r length(keep)` genes had enough counts to be tested (`edgeR::filterByExpr`).

# 4. Sample relationships

```{r transform}
# Variance-stabilising transformation (VST) makes counts comparable across genes for
# PCA and heatmaps. vst() needs >= 1000 genes; smaller data use the slower full version,
# and if that fails we fall back to log2(normalised counts + 1).
dds <- DESeqDataSetFromMatrix(counts_f, samples, design_formula)
dds <- estimateSizeFactors(dds)
transform_used <- "VST"
vsd_mat <- tryCatch(
  assay(if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)),
  error = function(e) {
    transform_used <<- "log2(normalised counts + 1)"
    log2(counts(dds, normalized = TRUE) + 1)
  })
```

Transformation used for this section: **`r transform_used`**.

```{r pca-function}
# PCA on the 500 most variable genes. Colour = condition, shape = first covariate (e.g. batch).
covariates <- setdiff(design_vars, "condition")
plot_pca <- function(mat, title) {
  top <- head(order(apply(mat, 1, var), decreasing = TRUE), min(500, nrow(mat)))
  pca <- prcomp(t(mat[top, , drop = FALSE]))
  var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  df <- cbind(as.data.frame(pca$x[, 1:2]), samples)
  p <- ggplot(df, aes(PC1, PC2, colour = condition, label = sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3, show.legend = FALSE) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2]))
  if (length(covariates) > 0 && is.factor(df[[covariates[1]]])) p <- p + aes(shape = .data[[covariates[1]]])
  list(plot = p, var_pct = var_pct)
}
pca_raw <- plot_pca(vsd_mat, "PCA")
pca_raw$plot
```

```{r scree, fig.height = 3}
# How much variation each principal component explains.
n_pc <- min(4, length(pca_raw$var_pct))
ggplot(data.frame(PC = paste0("PC", seq_len(n_pc)), pct = pca_raw$var_pct[seq_len(n_pc)]), aes(PC, pct)) +
  geom_col() + labs(x = NULL, y = "% variance explained")
```

```{r distance-heatmap, fig.height = 5}
# Euclidean distances between samples: replicates should cluster together.
pheatmap(as.matrix(dist(t(vsd_mat))), annotation_col = samples[, design_vars, drop = FALSE],
         main = "Sample-to-sample distances")
```

```{r batch-corrected-pca}
# If the design has covariates, show the PCA after removing them (for display only;
# DE uses the covariates in the model instead of corrected values).
if (length(covariates) > 0) {
  factor_covs <- covariates[vapply(covariates, function(v) is.factor(samples[[v]]), logical(1))]
  numeric_covs <- setdiff(covariates, factor_covs)
  corrected <- removeBatchEffect(
    vsd_mat,
    batch = if (length(factor_covs) >= 1) samples[[factor_covs[1]]],
    batch2 = if (length(factor_covs) >= 2) samples[[factor_covs[2]]],
    covariates = if (length(numeric_covs) > 0) as.matrix(samples[, numeric_covs, drop = FALSE]),
    design = model.matrix(~ condition, samples))
  plot_pca(corrected, paste("PCA after removing", paste(covariates, collapse = " + ")))$plot
}
```

# 5. Differential expression {.tabset}

```{r differential-expression}
# DESeq2 (Wald test, apeglm-shrunk fold changes) and edgeR (quasi-likelihood F-test)
# run on the same filtered counts. A gene is significant when
# padj < padj_cutoff and |log2FC| > lfc_cutoff (unshrunk log2FC, for both methods).
is_significant <- function(x) !is.na(x$padj) & x$padj < params$padj_cutoff & abs(x$log2FC) > params$lfc_cutoff
de <- list()
if (can_test) {
  dds <- DESeq(dds, quiet = TRUE)
  # edgeR: one column per condition (~ 0 + condition + covariates), so contrasts are differences of columns.
  mm_edger <- model.matrix(as.formula(paste(c("~ 0 + condition", covariates), collapse = " + ")), samples)
  y <- estimateDisp(calcNormFactors(DGEList(counts_f)), mm_edger)
  fit <- glmQLFit(y, mm_edger)

  for (cn in names(contrast_list)) {
    a <- contrast_list[[cn]][1]
    b <- contrast_list[[cn]][2]
    # DESeq2: make B the reference level so the coefficient is exactly "condition_A_vs_B".
    dds_b <- dds
    dds_b$condition <- relevel(dds_b$condition, ref = b)
    dds_b <- nbinomWaldTest(dds_b, quiet = TRUE)
    coef <- paste0("condition_", a, "_vs_", b)
    res <- results(dds_b, name = coef, alpha = params$padj_cutoff)
    shrunk <- lfcShrink(dds_b, coef = coef, res = res, type = "apeglm", quiet = TRUE)
    d <- data.frame(gene_id = rownames(res), gene_name = unname(symbol_of[rownames(res)]),
                    baseMean = res$baseMean, log2FC = res$log2FoldChange,
                    log2FC_shrunk = shrunk$log2FoldChange, pvalue = res$pvalue, padj = res$padj)
    # edgeR
    con <- makeContrasts(contrasts = sprintf("condition%s - condition%s", a, b), levels = mm_edger)
    tt <- topTags(glmQLFTest(fit, contrast = con), n = Inf, sort.by = "none")$table
    e <- data.frame(gene_id = rownames(tt), gene_name = unname(symbol_of[rownames(tt)]),
                    logCPM = tt$logCPM, log2FC = tt$logFC, pvalue = tt$PValue, padj = tt$FDR)
    d$significant <- is_significant(d)
    e$significant <- is_significant(e)
    de[[cn]] <- list(deseq2 = d[order(d$padj), ], edger = e[order(e$padj), ])
  }
}

# Save full result tables (one CSV per contrast and method).
dir.create(params$tables_dir, showWarnings = FALSE, recursive = TRUE)
for (cn in names(de)) {
  write.csv(de[[cn]]$deseq2, file.path(params$tables_dir, paste0(cn, "_deseq2.csv")), row.names = FALSE)
  write.csv(de[[cn]]$edger, file.path(params$tables_dir, paste0(cn, "_edger.csv")), row.names = FALSE)
}
```

```{r contrast-tabs, echo = FALSE, results = "asis"}
# One tab per contrast, built from the template in _contrast.Rmd.
if (can_test) {
  for (cn in names(de)) {
    cat(knitr::knit_child(file.path(REPORT_DIR, "_contrast.Rmd"), envir = environment(), quiet = TRUE), sep = "\n")
  }
} else {
  cat("Differential expression was skipped (see Warnings).\n")
}
```

# 6. DESeq2 vs edgeR

```{r concordance}
# How well the two methods agree: overlap of significant genes and fold-change correlation.
if (can_test) {
  conc <- do.call(rbind, lapply(names(de), function(cn) {
    m <- merge(de[[cn]]$deseq2[, c("gene_id", "log2FC", "significant")],
               de[[cn]]$edger[, c("gene_id", "log2FC", "significant")],
               by = "gene_id", suffixes = c("_deseq2", "_edger"))
    m$contrast <- cn
    m$category <- ifelse(m$significant_deseq2 & m$significant_edger, "both",
                  ifelse(m$significant_deseq2, "DESeq2 only",
                  ifelse(m$significant_edger, "edgeR only", "neither")))
    m
  }))
  concordance <- do.call(rbind, lapply(split(conc, conc$contrast), function(m) data.frame(
    contrast = m$contrast[1],
    both = sum(m$category == "both"),
    deseq2_only = sum(m$category == "DESeq2 only"),
    edger_only = sum(m$category == "edgeR only"),
    spearman_log2FC = round(cor(m$log2FC_deseq2, m$log2FC_edger, method = "spearman", use = "complete.obs"), 3))))
  knitr::kable(concordance, row.names = FALSE)
}
```

```{r concordance-plot, fig.width = 8}
# Each point is a gene: DESeq2 fold change (x) against edgeR fold change (y).
if (can_test) {
  ggplot(conc, aes(log2FC_deseq2, log2FC_edger, colour = category)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_abline(linetype = "dashed") +
    scale_colour_manual(values = c(both = "firebrick", `DESeq2 only` = "darkorange",
                                   `edgeR only` = "steelblue", neither = "grey70")) +
    facet_wrap(~ contrast) +
    labs(x = "DESeq2 log2FC", y = "edgeR log2FC", colour = NULL)
}
```

# 7. GO enrichment

```{r go-enrichment, results = "asis", fig.width = 8, fig.height = 6}
# Over-representation of GO terms (BP, MF, CC) among DESeq2 significant genes,
# up- and down-regulated separately. Background = all genes that passed filtering.
orgdb <- switch(params$genome, GRCh38 = "org.Hs.eg.db", GRCm39 = "org.Mm.eg.db", NULL)
go_enabled <- can_test && !isTRUE(params$skip_enrichment) && !is.null(orgdb)
run_go <- function(genes) {
  if (length(genes) < MIN_GO_GENES) return(NULL)
  ego <- tryCatch(
    clusterProfiler::enrichGO(genes, OrgDb = orgdb, keyType = "ENSEMBL", ont = "ALL",
                              universe = rownames(counts_f), pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) NULL)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) NULL else ego
}
if (!go_enabled) {
  cat(if (isTRUE(params$skip_enrichment)) "GO enrichment was skipped (--skip_enrichment).\n"
      else if (is.null(orgdb)) "GO enrichment is only available for GRCh38 and GRCm39.\n"
      else "Differential expression was skipped, so there is nothing to test.\n")
} else {
  for (cn in names(de)) {
    d <- de[[cn]]$deseq2
    for (direction in c("up", "down")) {
      genes <- d$gene_id[d$significant & (if (direction == "up") d$log2FC > 0 else d$log2FC < 0)]
      ego <- run_go(genes)
      cat("\n\n## ", cn, ": ", direction, "-regulated (", length(genes), " genes)\n\n", sep = "")
      if (is.null(ego)) {
        cat("No enriched GO terms (or fewer than ", MIN_GO_GENES, " significant genes).\n\n", sep = "")
      } else {
        write.csv(as.data.frame(ego), file.path(params$tables_dir, paste0(cn, "_go_", direction, ".csv")), row.names = FALSE)
        print(enrichplot::dotplot(ego, showCategory = 10, split = "ONTOLOGY") +
                facet_grid(ONTOLOGY ~ ., scales = "free_y") + ggtitle(paste(cn, direction)))
        cat("\n\n")
      }
    }
  }
}
```

# 8. Exported files

```{r exports, results = "asis"}
# Every CSV written by this report (in the tables/ folder next to the report).
files <- list.files(params$tables_dir)
cat(if (length(files) > 0) paste0("- `tables/", files, "`", collapse = "\n") else "No tables were written.", "\n")
```

## R session

```{r session-info}
# Package versions, for reproducibility.
sessionInfo()
```
````

- [ ] **Step 4: Write `report/_contrast.Rmd`**

````markdown
## `r cn`

```{r}
# Number of significant genes per method and direction.
r <- de[[cn]]
n_sig <- sapply(r[c("deseq2", "edger")], function(x) c(
  up = sum(x$significant & x$log2FC > 0), down = sum(x$significant & x$log2FC < 0)))
knitr::kable(n_sig, col.names = c("DESeq2", "edgeR"))
```

Significant = padj < `r params$padj_cutoff` and |log2FC| > `r params$lfc_cutoff`. Positive log2FC = higher in **`r contrast_list[[cn]][1]`** than in **`r contrast_list[[cn]][2]`**.

```{r, fig.width = 10, fig.height = 4}
# MA plot: mean expression (x) against fold change (y). Red = significant.
# DESeq2 shows apeglm-shrunk fold changes, which pull noisy low-count genes towards 0.
ma <- rbind(
  data.frame(method = "DESeq2 (shrunk log2FC)", mean_expr = log10(r$deseq2$baseMean + 1),
             log2FC = r$deseq2$log2FC_shrunk, significant = r$deseq2$significant),
  data.frame(method = "edgeR", mean_expr = r$edger$logCPM,
             log2FC = r$edger$log2FC, significant = r$edger$significant))
ggplot(ma, aes(mean_expr, log2FC, colour = significant)) +
  geom_point(size = 0.8, alpha = 0.6) + geom_hline(yintercept = 0) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
  facet_wrap(~ method, scales = "free_x") +
  labs(x = "Mean expression (DESeq2: log10 baseMean, edgeR: logCPM)", y = "log2 fold change")
```

```{r, fig.width = 10, fig.height = 4.5}
# Volcano plot: fold change (x) against significance (y). Top 10 genes per method are labelled.
cols <- c("gene_name", "log2FC", "pvalue", "padj", "significant")
volcano <- rbind(cbind(method = "DESeq2", r$deseq2[, cols]), cbind(method = "edgeR", r$edger[, cols]))
volcano <- volcano[!is.na(volcano$pvalue), ]
labels <- do.call(rbind, lapply(split(volcano, volcano$method), function(x) {
  s <- x[x$significant, ]
  head(s[order(s$padj), ], 10)
}))
ggplot(volcano, aes(log2FC, -log10(pvalue), colour = significant)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1) * params$lfc_cutoff, linetype = "dashed") +
  geom_text(data = labels, aes(label = gene_name), size = 3, vjust = -0.6, colour = "black", check_overlap = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
  facet_wrap(~ method)
```

```{r, fig.height = 7}
# Heatmap of up to 50 top DESeq2 genes (row z-scores of transformed counts).
top <- head(r$deseq2$gene_id[r$deseq2$significant], 50)
if (length(top) > 0) top <- top[apply(vsd_mat[top, , drop = FALSE], 1, sd) > 0]  # constant rows break scaling
if (length(top) >= 2) {
  z <- t(scale(t(vsd_mat[top, , drop = FALSE])))
  rownames(z) <- make.unique(unname(symbol_of[top]))
  pheatmap(z, annotation_col = samples[, design_vars, drop = FALSE], fontsize_row = 7,
           main = paste("Top significant genes:", cn))
} else {
  cat("Fewer than 2 significant genes, so no heatmap.\n")
}
```

```{r}
# DESeq2 results (top 1000 by padj; the full table is in tables/).
DT::formatSignif(DT::datatable(head(r$deseq2, 1000), rownames = FALSE, filter = "top",
                               caption = "DESeq2 results"),
                 c("baseMean", "log2FC", "log2FC_shrunk", "pvalue", "padj"), 3)
```

```{r}
# edgeR results (top 1000 by padj; the full table is in tables/).
DT::formatSignif(DT::datatable(head(r$edger, 1000), rownames = FALSE, filter = "top",
                               caption = "edgeR results"),
                 c("logCPM", "log2FC", "pvalue", "padj"), 3)
```
````

- [ ] **Step 5: Standalone render against the Task 5 results**

Run:
```bash
mkdir -p results/standalone
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/w -w /w/report bulk-rnaseq-lite-report:1.0 \
  Rscript -e 'rmarkdown::render("rnaseq_report.Rmd", output_dir = "../results/standalone", params = list(design = "~ batch + condition", tables_dir = "../results/standalone/tables"))'
ls results/standalone/
```
Expected: `rnaseq_report.html` and `tables/` containing `treated_vs_control_deseq2.csv` and `treated_vs_control_edger.csv`. Open the HTML and check that every section renders and the Run summary shows `reverse (auto-detected)`.

- [ ] **Step 6: Review Focus 3 — nothing significant**

Run:
```bash
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/w -w /w/report bulk-rnaseq-lite-report:1.0 \
  Rscript -e 'rmarkdown::render("rnaseq_report.Rmd", output_file = "none_sig.html", output_dir = "../results/standalone", params = list(lfc_cutoff = 100, tables_dir = "../results/standalone/tables_none"))'
grep -c "Fewer than 2 significant genes" results/standalone/none_sig.html
```
Expected: renders with no error; count ≥ 1.

- [ ] **Step 7: Review Focus 4 — batch confounded with condition**

Run:
```bash
awk -F, 'BEGIN{OFS=","} NR==1{print; next} {$5=$4; print}' results/pipeline_info/samplesheet.csv > results/standalone/confounded.csv
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/w -w /w/report bulk-rnaseq-lite-report:1.0 \
  Rscript -e 'rmarkdown::render("rnaseq_report.Rmd", output_file = "confounded.html", output_dir = "../results/standalone", params = list(samplesheet = "../results/standalone/confounded.csv", design = "~ batch + condition", tables_dir = "../results/standalone/tables_conf"))'
grep -c "confounded with condition" results/standalone/confounded.html
```
Expected: renders; count ≥ 1.

- [ ] **Step 8: Write `modules/report.nf`**

```groovy
// Render the R Markdown report into one self-contained HTML file.
// The same Rmd can be knitted by hand; see docs/usage.md.
process REPORT {
    label 'process_medium'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path report_dir, stageAs: 'report_src'
    path counts
    path samplesheet, stageAs: 'input_samplesheet.csv'
    path gene_names
    path qc_summary
    path strand_report
    path versions
    val strandedness

    output:
    path 'report.html', emit: html
    path 'tables', emit: tables

    script:
    def release = params.gtf ? 'custom reference' : params.ensembl_release
    def contrasts = params.contrasts ?: ''
    def skip_go = params.skip_enrichment ? 'TRUE' : 'FALSE'
    """
    export HOME=\$PWD
    cp -rL report_src report_work
    Rscript -e "rmarkdown::render('report_work/rnaseq_report.Rmd', \\
        output_file = 'report.html', output_dir = getwd(), \\
        knit_root_dir = getwd(), intermediates_dir = getwd(), \\
        params = list(counts = '${counts}', samplesheet = '${samplesheet}', \\
                      gene_names = '${gene_names}', qc_summary = '${qc_summary}', \\
                      strand_report = '${strand_report}', versions = '${versions}', \\
                      multiqc_link = '../multiqc/multiqc_report.html', \\
                      design = '${params.design}', contrasts = '${contrasts}', \\
                      padj_cutoff = ${params.padj_cutoff}, lfc_cutoff = ${params.lfc_cutoff}, \\
                      genome = '${params.genome}', ensembl_release = '${release}', \\
                      strandedness = '${strandedness}', skip_enrichment = ${skip_go}, \\
                      tables_dir = 'tables'))"
    """
}
```

- [ ] **Step 9: Wire REPORT into `main.nf`**

Add to the includes:

```groovy
include { REPORT        } from './modules/report'
```

Append at the end of the `workflow` block:

```groovy
    // ---- Downstream analysis and HTML report -----------------------------
    REPORT(
        file("${projectDir}/report"),
        MERGE_COUNTS.out,
        samplesheet_copy,
        GENE_NAMES.out,
        QC_SUMMARY.out,
        strand_report,
        versions,
        strandedness
    )
```

- [ ] **Step 10: Lint, run the full test profile and the end-to-end check**

Run: `nextflow lint main.nf modules/ && nextflow run main.nf -profile test,docker -resume && python3 tests/check_results.py`
Expected: pipeline completes; `ok deseq2: ...`, `ok edger: ...`, `ok report.html (... KB), strandedness=reverse`.

- [ ] **Step 11: Commit**

```bash
git add report/ modules/report.nf main.nf tests/check_results.py
git commit -m "feat: add R Markdown report with DESeq2, edgeR, GO and QC sections

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Documentation and example params file

**Files:**
- Create: `README.md`, `docs/usage.md`, `docs/output.md`, `params.yaml`

**Interfaces:**
- Consumes: parameter names and defaults from `nextflow.config` (Task 4), output paths (Tasks 5 and 7), Rmd params (Task 7).

- [ ] **Step 1: Write `params.yaml`**

```yaml
# Example run configuration. Copy it, edit it, then run:
#   nextflow run main.nf -profile docker -params-file params.yaml
# Anything given on the command line (--name value) overrides this file.

input: samplesheet.csv       # required: see docs/usage.md#samplesheet
outdir: results

genome: GRCh38               # GRCh38 (human) or GRCm39 (mouse)
ensembl_release: 112

strandedness: auto           # auto | forward | reverse | unstranded
read_length: 100             # used for the STAR index; 100 works for most data

design: "~ condition"        # e.g. "~ batch + condition"
contrasts: null              # e.g. "treated_vs_control,KO_vs_WT"; null = every condition vs the first
padj_cutoff: 0.05
lfc_cutoff: 0.58             # log2 fold change; 0.58 ~ 1.5-fold

skip_trimming: false
skip_enrichment: false
save_bam: false
```

- [ ] **Step 2: Write `README.md`**

````markdown
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
````

- [ ] **Step 3: Write `docs/usage.md`**

````markdown
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
````

- [ ] **Step 4: Write `docs/output.md`**

````markdown
# Outputs

```
results/
├── fastqc/               FastQC reports for raw reads (<sample>_fastqc.html)
├── fastp/                Trimming reports (<sample>.fastp.html / .json)
├── star/                 <sample>.Log.final.out (mapping stats), <sample>.ReadsPerGene.out.tab,
│                         strandedness.txt + strandedness_mqc.tsv (auto mode), BAMs with --save_bam
├── featurecounts/        counts.tsv (genes x samples, raw counts), gene_names.tsv,
│   └── per_sample/       featureCounts output and summary per sample
├── qc/qc_summary.tsv     Reads, % retained, % uniquely mapped, % assigned per sample
├── multiqc/              multiqc_report.html: all upstream QC in one page
├── report/               report.html + tables/ (CSV per contrast)
├── reference/            STAR index, only when --fasta/--gtf were given
└── pipeline_info/        samplesheet.csv, versions.tsv, timeline/trace/execution reports
```

`counts.tsv` holds **raw** counts: use it for DESeq2/edgeR, not for plotting expression directly.

## Result tables (`report/tables/`)

`<A>_vs_<B>_deseq2.csv`: `gene_id, gene_name, baseMean, log2FC, log2FC_shrunk, pvalue, padj, significant`
`<A>_vs_<B>_edger.csv`: `gene_id, gene_name, logCPM, log2FC, pvalue, padj, significant`
`<A>_vs_<B>_go_up.csv` / `_go_down.csv`: enriched GO terms (only written when terms are found)

Positive log2FC = higher in A. `significant` = padj < `padj_cutoff` and |log2FC| > `lfc_cutoff`
(unshrunk log2FC for both methods). `log2FC_shrunk` (apeglm) is better for ranking genes and
for plots because noisy low-count genes are pulled towards 0.

## Reading the report

**Warnings and checks.** Read first. Lists anything that changed the analysis: a design that couldn't be
estimated (replaced by `~ condition`), groups with < 2 replicates, low gene assignment.

**1. Run summary.** Genome, Ensembl release, strandedness used, design, contrasts, cutoffs, tool versions.

**2. Upstream QC.** Per-sample table and bars. Uniquely mapped is typically 70–90% for good human/mouse
data; below 60% is flagged. Assigned to genes below 50% is flagged, usually a strandedness or annotation issue.

**3. Count QC.** Library sizes should be of the same order across samples. The log-CPM boxes
should look similar; one sample far off is worth checking.

**4. Sample relationships.**
- *PCA*: each dot is a sample. Replicates of a condition should sit together; if samples split by batch
  instead, batch is a strong effect and belongs in `--design`.
- *Scree*: how much variation each PC captures.
- *Distance heatmap*: darker = more similar. Look for outliers.
- *PCA after removing covariates*: shows what the data look like once batch is accounted for. For display only.

**5. Differential expression** (one tab per contrast).
- *MA plot*: expression (x) vs fold change (y). Red points are significant.
- *Volcano*: fold change (x) vs −log10 p-value (y). Top genes labelled.
- *Heatmap*: top significant genes, scaled per gene (z-score); columns should group by condition.
- *Tables*: searchable, top 1000 genes; full tables in `tables/`.

**6. DESeq2 vs edgeR.** Genes significant in both methods are the most robust. Fold changes
normally correlate strongly (Spearman > 0.9).

**7. GO enrichment.** Over-represented Biological Process / Molecular Function / Cellular Component
terms among significant DESeq2 genes, up and down separately. Dot size = number of genes,
colour = adjusted p-value. Fewer than 5 significant genes: skipped.
````

- [ ] **Step 5: Check every parameter in `nextflow.config` is documented**

Run: `for p in $(sed -n '/^params {/,/^}/p' nextflow.config | grep -oP '^\s{4}\K[a-z_]+(?=\s+=)' | grep -v containers | sort -u); do grep -q -- "--$p" docs/usage.md || echo "undocumented: $p"; done`
Expected: `undocumented: low_memory` only (it's set through the profile, and the Profiles section covers it). Nothing else.

- [ ] **Step 6: Commit**

```bash
git add README.md docs/usage.md docs/output.md params.yaml
git commit -m "docs: add README, usage, output guide and example params file

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Continuous integration

**Files:**
- Create: `.github/workflows/ci.yml`

**Interfaces:**
- Consumes: `tests/test_bin.py`, `tests/validation.sh`, `tests/check_results.py`, `docker/Dockerfile`, test profile.

- [ ] **Step 1: Write `.github/workflows/ci.yml`**

```yaml
name: CI
on:
  push:
  pull_request:

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-java@v4
        with:
          distribution: temurin
          java-version: "21"
      - uses: nf-core/setup-nextflow@v2
      - name: Python helper tests
        run: python3 tests/test_bin.py
      - name: Lint Nextflow code
        run: nextflow lint main.nf modules/ nextflow.config
      - name: Input validation tests
        run: bash tests/validation.sh
      - name: Build report image
        run: docker build -t bulk-rnaseq-lite-report:1.0 docker/
      - name: Run pipeline on test data
        run: nextflow run main.nf -profile test,docker
      - name: Check results against simulated truth
        run: python3 tests/check_results.py
      - name: Upload report
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: test-report
          path: |
            results/report/
            results/multiqc/multiqc_report.html
            .nextflow.log
```

- [ ] **Step 2: Commit and push**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: run helper, validation and end-to-end tests on push

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push origin main
```

- [ ] **Step 3: Watch the run**

Run: `gh run watch --exit-status $(gh run list --limit 1 --json databaseId -q '.[0].databaseId')`
Expected: all steps green. If a step fails, read the log with `gh run view --log-failed` and fix the cause before continuing.

- [ ] **Step 4: Add the CI badge to the README**

Insert under the title line in `README.md`:

```markdown
[![CI](https://github.com/DimAnge/bulk-rnaseq-lite-nf/actions/workflows/ci.yml/badge.svg)](https://github.com/DimAnge/bulk-rnaseq-lite-nf/actions/workflows/ci.yml)
```

(The badge only renders for viewers with access while the repo is private.)

```bash
git add README.md
git commit -m "docs: add CI badge

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push origin main
```

---

### Task 10: Real-data smoke test (user, on the HPC server)

This can't run on the dev laptop (15 GB RAM) with a full index. It checks the parts no automated test covers: Ensembl download, full-size index build and real-data runtime.

- [ ] **Step 1: On the HPC server, clone and build the image**

```bash
git clone git@github.com:DimAnge/bulk-rnaseq-lite-nf.git && cd bulk-rnaseq-lite-nf
docker build -t bulk-rnaseq-lite-report:1.0 docker/
```

- [ ] **Step 2: Run on a small real dataset (4–6 samples, human or mouse)**

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --genome GRCm39 --max_cpus 16 --max_memory 64.GB
```

Expected: `DOWNLOAD_REFS` and `STAR_INDEX` run once, `~/.rnaseq-refs/GRCm39/ensembl_112/` holds `genome.fa`, `genes.gtf`, `star_index/`; the report opens and the strandedness call matches the library kit.

- [ ] **Step 3: Run again and confirm no re-download or re-index**

Run the same command with `--outdir results2`. Expected: neither `DOWNLOAD_REFS` nor `STAR_INDEX` executes.

- [ ] **Step 4: Report anything odd back into the repo as issues or fixes before making it public.**
