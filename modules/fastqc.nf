// Quality report for the raw reads.
// Inputs are staged into input1/, input2/ so two samples with the same FASTQ file name
// never collide; the links give each report a sample-based name.
process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads, stageAs: 'input?/*')

    output:
    path '*_fastqc.zip', emit: zip
    path '*_fastqc.html', emit: html

    script:
    def r = [reads].flatten()   // always a list, for single- and paired-end
    def links = meta.single_end
        ? "ln -s ${r[0]} ${meta.id}.fastq.gz"
        : "ln -s ${r[0]} ${meta.id}_1.fastq.gz && ln -s ${r[1]} ${meta.id}_2.fastq.gz"
    """
    ${links}
    fastqc --threads ${task.cpus} --quiet ${meta.id}*.fastq.gz
    """
}
