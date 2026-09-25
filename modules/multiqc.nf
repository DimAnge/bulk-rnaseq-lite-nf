// One HTML page with FastQC, fastp, STAR, featureCounts and strandedness results.
process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc . --force --filename multiqc_report.html
    """
}
