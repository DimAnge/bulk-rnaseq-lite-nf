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
