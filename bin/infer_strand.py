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
