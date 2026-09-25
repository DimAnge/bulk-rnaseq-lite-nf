// Download the Ensembl genome FASTA and GTF once. storeDir makes Nextflow skip this
// process entirely when both files already exist in the cache folder.
process DOWNLOAD_REFS {
    tag "${store_dir}"
    storeDir "${store_dir}"

    input:
    val fasta_url
    val gtf_url
    val store_dir

    output:
    path 'genome.fa', emit: fasta
    path 'genes.gtf', emit: gtf

    script:
    """
    python3 -c "import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], 'genome.fa.gz')" ${fasta_url}
    python3 -c "import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], 'genes.gtf.gz')" ${gtf_url}
    gunzip genome.fa.gz genes.gtf.gz
    """
}
