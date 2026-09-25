// Adapter and quality trimming. fastp detects adapters by itself and writes a JSON
// report that MultiQC and the QC summary read.
process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(meta), path(reads, stageAs: 'input?/*')

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html", emit: html

    script:
    def r = [reads].flatten()
    def io = meta.single_end
        ? "--in1 ${r[0]} --out1 ${meta.id}.trim.fastq.gz"
        : "--in1 ${r[0]} --in2 ${r[1]} --out1 ${meta.id}_1.trim.fastq.gz --out2 ${meta.id}_2.trim.fastq.gz --detect_adapter_for_pe"
    """
    fastp ${io} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html
    """
}
