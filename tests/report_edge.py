#!/usr/bin/env python3
"""Render the report standalone on edge-case inputs built from a finished test run.

Needs: results/ from `nextflow run main.nf -profile test,docker` and the report image.
Usage: python3 tests/report_edge.py [results_dir]
"""
import csv
import html
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "results"
IMAGE = "bulk-rnaseq-lite-report:1.0"


def read_tsv(path):
    with open(path) as fh:
        return list(csv.reader(fh, delimiter="\t"))


def write_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)


def render(tmp, name, **params):
    """Render the Rmd with params (paths relative to /data = tmp); return visible text (code blocks removed)."""
    base = {"counts": "/r/featurecounts/counts.tsv", "samplesheet": "/r/pipeline_info/samplesheet.csv",
            "gene_names": "/r/featurecounts/gene_names.tsv", "qc_summary": "/r/qc/qc_summary.tsv",
            "strand_report": "/r/star/strandedness_mqc.tsv", "versions": "/r/pipeline_info/versions.tsv",
            "design": "~ batch + condition", "tables_dir": f"/data/{name}_tables"}
    base.update(params)
    r_params = ", ".join(f"{k} = {v}" if isinstance(v, (int, float)) else f"{k} = '{v}'" for k, v in base.items())
    cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
           "-v", f"{ROOT / 'report'}:/report:ro", "-v", f"{RESULTS}:/r:ro", "-v", f"{tmp}:/data", "-w", "/data", IMAGE,
           "Rscript", "-e", f"file.copy(list.files('/report', full.names = TRUE), '/data'); "
           f"rmarkdown::render('/data/rnaseq_report.Rmd', output_file = '{name}.html', output_dir = '/data', "
           f"quiet = TRUE, params = list({r_params}))"]
    out = subprocess.run(cmd, capture_output=True, text=True)
    assert out.returncode == 0, f"{name}: render failed\n{out.stderr[-2000:]}"
    page = (Path(tmp) / f"{name}.html").read_text()
    page = re.sub(r'<pre class="r">.*?</pre>|<script.*?</script>|<style.*?</style>', "", page, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", page)))


def expect(text, pattern, name):
    assert re.search(pattern, text), f"{name}: expected /{pattern}/ in the rendered report"
    print(f"ok {name}: /{pattern}/")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        counts = read_tsv(RESULTS / "featurecounts/counts.tsv")
        names = read_tsv(RESULTS / "featurecounts/gene_names.tsv")
        sheet = list(csv.reader(open(RESULTS / "pipeline_info/samplesheet.csv")))

        # 1. Nothing significant.
        t = render(tmp, "none_sig", lfc_cutoff=100)
        expect(t, r"Fewer than 2 significant genes, so no heatmap", "none_sig")

        # 2. Batch confounded with condition -> falls back to ~ condition with a warning.
        conf = [sheet[0]] + [r[:4] + [r[3]] for r in sheet[1:]]
        csv.writer(open(f"{tmp}/confounded.csv", "w")).writerows(conf)
        t = render(tmp, "confounded", samplesheet="/data/confounded.csv")
        expect(t, r"Warning: The design .~ batch \+ condition. cannot be estimated", "confounded")

        # 3. Numeric-looking sample IDs (01..06) and batch coded 1/2/3.
        ids = {old: f"{i + 1:02d}" for i, old in enumerate(counts[0][1:])}
        write_tsv(f"{tmp}/num_counts.tsv", [["gene_id"] + [ids[c] for c in counts[0][1:]]] + counts[1:])
        num = [sheet[0]] + [[ids[r[0]]] + r[1:4] + [str(i % 3 + 1)] for i, r in enumerate(sheet[1:])]
        csv.writer(open(f"{tmp}/num_sheet.csv", "w")).writerows(num)
        t = render(tmp, "numeric", counts="/data/num_counts.tsv", samplesheet="/data/num_sheet.csv")
        expect(t, r"Column .batch. has only 3 whole-number values", "numeric")  # pandoc curls the quotes
        expect(t, r"treated_vs_control", "numeric")

        # 4. GENCODE-style versioned gene IDs still get GO-annotated.
        write_tsv(f"{tmp}/ver_counts.tsv", [counts[0]] + [[r[0] + ".7"] + r[1:] for r in counts[1:]])
        write_tsv(f"{tmp}/ver_names.tsv", [names[0]] + [[r[0] + ".7"] + r[1:] for r in names[1:]])
        t = render(tmp, "versioned", counts="/data/ver_counts.tsv", gene_names="/data/ver_names.tsv")
        expect(t, r"[1-9][0-9]* of [0-9]+ tested genes were found in org\.Hs\.eg\.db", "versioned")

        # 5. Human data with --genome GRCm39 -> clear mismatch message instead of "no terms".
        t = render(tmp, "wrong_species", genome="GRCm39")
        expect(t, r"Gene IDs do not match org\.Mm\.eg\.db", "wrong_species")
    print("all report edge cases passed")


if __name__ == "__main__":
    main()
