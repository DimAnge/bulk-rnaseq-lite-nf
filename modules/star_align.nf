// Align reads with STAR. --quantMode GeneCounts also writes ReadsPerGene.out.tab,
// which INFER_STRAND uses to detect the library strandedness.
process STAR_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/star", mode: 'copy',
        saveAs: { fn -> (fn.endsWith('.bam') && !params.save_bam) ? null : fn }

    input:
    tuple val(meta), path(reads, stageAs: 'input?/*')
    path index

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${meta.id}.Log.final.out", emit: log
    path "${meta.id}.ReadsPerGene.out.tab", emit: gene_counts

    script:
    def files = [reads].flatten().join(' ')
    """
    STAR --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${files} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS nM \\
        --quantMode GeneCounts \\
        --limitBAMsortRAM ${task.memory.toBytes().intdiv(2)}
    """
}
