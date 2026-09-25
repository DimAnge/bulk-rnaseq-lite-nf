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
