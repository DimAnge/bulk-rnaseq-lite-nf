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
