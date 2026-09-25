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
