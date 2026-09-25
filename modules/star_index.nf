// Build the STAR genome index once and keep it in the cache (storeDir).
// -profile low_memory builds a sparse index that needs ~16 GB instead of ~32 GB.
process STAR_INDEX {
    tag "${params.low_memory ? 'sparse' : 'full'} index"
    storeDir "${store_dir}"

    input:
    path fasta
    path gtf
    val store_dir

    output:
    path "${params.low_memory ? 'star_index_sparse3' : 'star_index'}", emit: index

    script:
    def index_dir = params.low_memory ? 'star_index_sparse3' : 'star_index'
    def sparse = params.low_memory ? '--genomeSAsparseD 3' : ''
    """
    # STAR cannot read gzipped references
    if [[ ${fasta} == *.gz ]]; then gunzip -c ${fasta} > ref.fa; else ln -s ${fasta} ref.fa; fi
    if [[ ${gtf} == *.gz ]]; then gunzip -c ${gtf} > ref.gtf; else ln -s ${gtf} ref.gtf; fi

    # Small genomes need a smaller suffix-array index: min(14, log2(genome length)/2 - 1)
    genome_len=\$(grep -v '>' ref.fa | tr -d '\\n' | wc -c)
    sa_bases=\$(awk -v len="\$genome_len" 'BEGIN { n = int(log(len) / log(2) / 2 - 1); print (n < 14 ? n : 14) }')

    mkdir ${index_dir}
    STAR --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index_dir} \\
        --genomeFastaFiles ref.fa \\
        --sjdbGTFfile ref.gtf \\
        --sjdbOverhang ${params.read_length - 1} \\
        --genomeSAindexNbases \$sa_bases \\
        --limitGenomeGenerateRAM ${task.memory.toBytes()} \\
        ${sparse}
    """
}
