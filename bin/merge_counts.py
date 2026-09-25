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
