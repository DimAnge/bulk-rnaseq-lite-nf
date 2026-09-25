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
