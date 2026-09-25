// Count reads per gene (exons summarised by gene_id). Runs per sample because one
// featureCounts call cannot mix paired-end (-p) and single-end BAMs.
process FEATURECOUNTS {
    tag "${meta.id}"
    publishDir "${params.outdir}/featurecounts/per_sample", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path gtf
    val strandedness

    output:
    path "${meta.id}.featureCounts.txt", emit: counts
    path "${meta.id}.featureCounts.txt.summary", emit: summary

    script:
    def strand_flag = [unstranded: 0, forward: 1, reverse: 2][strandedness]
    def paired = meta.single_end ? '' : '-p --countReadPairs'
    """
    featureCounts -T ${task.cpus} \\
        -a ${gtf} -t exon -g gene_id \\
        -s ${strand_flag} ${paired} \\
        -o ${meta.id}.featureCounts.txt \\
        ${bam}
    """
}
