// Render the R Markdown report into one self-contained HTML file.
// The same Rmd can be knitted by hand; see docs/usage.md.
process REPORT {
    label 'process_medium'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path report_dir, stageAs: 'report_src'
    path counts
    path samplesheet, stageAs: 'input_samplesheet.csv'
    path gene_names
    path qc_summary
    path strand_report
    path versions
    val strandedness

    output:
    path 'report.html', emit: html
    path 'tables', emit: tables

    script:
    def release = params.gtf ? 'custom reference' : params.ensembl_release
    def contrasts = params.contrasts ?: ''
    def skip_go = params.skip_enrichment ? 'TRUE' : 'FALSE'
    """
    export HOME=\$PWD
    cp -rL report_src report_work
    Rscript -e "rmarkdown::render('report_work/rnaseq_report.Rmd', \\
        output_file = 'report.html', output_dir = getwd(), \\
        knit_root_dir = getwd(), intermediates_dir = getwd(), \\
        params = list(counts = '${counts}', samplesheet = '${samplesheet}', \\
                      gene_names = '${gene_names}', qc_summary = '${qc_summary}', \\
                      strand_report = '${strand_report}', versions = '${versions}', \\
                      multiqc_link = '../multiqc/multiqc_report.html', \\
                      design = '${params.design}', contrasts = '${contrasts}', \\
                      padj_cutoff = ${params.padj_cutoff}, lfc_cutoff = ${params.lfc_cutoff}, \\
                      genome = '${params.genome}', ensembl_release = '${release}', \\
                      strandedness = '${strandedness}', skip_enrichment = ${skip_go}, \\
                      tables_dir = 'tables'))"
    """
}
