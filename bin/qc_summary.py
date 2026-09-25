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
