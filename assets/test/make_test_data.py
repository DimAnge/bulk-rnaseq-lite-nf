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
