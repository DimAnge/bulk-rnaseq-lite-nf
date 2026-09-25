#!/usr/bin/env nextflow
/*
 * bulk-rnaseq-lite-nf
 * FASTQ -> FastQC + fastp -> STAR -> featureCounts -> DESeq2 + edgeR + GO -> HTML report
 *
 * Read this file top to bottom: input checks first, then the workflow that wires the
 * modules together. Each module in modules/ is one tool.
 */

include { DOWNLOAD_REFS } from './modules/download_refs'
include { STAR_INDEX    } from './modules/star_index'
include { FASTQC        } from './modules/fastqc'
include { FASTP         } from './modules/fastp'
include { STAR_ALIGN    } from './modules/star_align'
include { INFER_STRAND  } from './modules/infer_strand'
include { FEATURECOUNTS } from './modules/featurecounts'
include { MERGE_COUNTS; GENE_NAMES; QC_SUMMARY } from './modules/helpers'
include { MULTIQC       } from './modules/multiqc'

// ---------------------------------------------------------------------------
// Input checks. They run before any process starts, so mistakes fail in seconds.
// ---------------------------------------------------------------------------

def validStrandedness() {
    return ['auto', 'forward', 'reverse', 'unstranded']
}

def supportedGenomes() {
    return ['GRCh38', 'GRCm39']
}

def checkParams() {
    if (!params.input) {
        error("Please provide a samplesheet with --input samplesheet.csv")
    }
    if (!(params.strandedness in validStrandedness())) {
        error("--strandedness must be one of ${validStrandedness().join(', ')} (got '${params.strandedness}')")
    }
    if (params.fasta && !params.gtf) {
        error("--fasta and --gtf must be given together")
    }
    if (params.gtf && !params.fasta && !params.star_index) {
        error("--gtf needs --fasta (to build an index) or --star_index")
    }
    if (!params.gtf && !(params.genome in supportedGenomes())) {
        error("--genome must be one of ${supportedGenomes().join(', ')} (got '${params.genome}'). For other species give --fasta and --gtf.")
    }
    if (!(params.padj_cutoff > 0 && params.padj_cutoff <= 1)) {
        error("--padj_cutoff must be between 0 and 1 (got ${params.padj_cutoff})")
    }
    if (params.lfc_cutoff < 0) {
        error("--lfc_cutoff must be 0 or more (got ${params.lfc_cutoff})")
    }
}

// Relative FASTQ paths are read relative to the samplesheet's folder.
def resolvePath(String p, base) {
    return (p.startsWith('/') || p.contains('://')) ? file(p) : base.resolve(p)
}

def readSamplesheet(String path) {
    def sheet = file(path)
    if (!sheet.exists()) {
        error("Samplesheet not found: ${path}")
    }
    // Trim keys and values; drop the invisible BOM that Excel adds to the first header.
    def rows = sheet.splitCsv(header: true).collect { row ->
        row.collectEntries { k, v -> [k.replace('﻿', '').trim(), (v ?: '').trim()] }
    }
    if (!rows) {
        error("Samplesheet ${path} has no samples")
    }
    def missing = ['sample', 'fastq_1', 'condition'].findAll { c -> !rows[0].containsKey(c) }
    if (missing) {
        error("Samplesheet is missing required column(s): ${missing.join(', ')}")
    }
    def names = rows.collect { r -> r.sample }
    def dups = names.findAll { n -> names.count(n) > 1 }.unique()
    if (dups) {
        error("Duplicate sample name(s) in samplesheet: ${dups.join(', ')}")
    }
    rows.each { r ->
        if (!r.sample.matches('[A-Za-z0-9_.\\-]+')) {
            error("Sample name '${r.sample}' may only contain letters, digits, '_', '.' and '-'")
        }
        if (!r.condition.matches('[A-Za-z][A-Za-z0-9_.]*') || r.condition.contains('_vs_')) {
            error("condition values must start with a letter, contain only letters, digits, '_' or '.', and not contain '_vs_' (got '${r.condition}' for sample ${r.sample})")
        }
        [r.fastq_1, r.fastq_2].findAll { p -> p }.each { p ->
            if (!(p.endsWith('.fastq.gz') || p.endsWith('.fq.gz'))) {
                error("FASTQ files must be gzipped and end in .fastq.gz or .fq.gz: ${p}")
            }
            if (!resolvePath(p, sheet.parent).exists()) {
                error("FASTQ file not found: ${p} (sample ${r.sample})")
            }
        }
    }
    return rows
}

def checkDesign(List rows) {
    def design = params.design.toString().trim()
    if (!design.startsWith('~')) {
        error("--design must start with '~', e.g. '~ batch + condition' (got '${design}')")
    }
    if (design.contains('*') || design.contains(':')) {
        error("--design interaction terms (* or :) are not supported (got '${design}')")
    }
    def vars = design.substring(1).tokenize('+').collect { v -> v.trim() }.findAll { v -> v }
    if (!vars || vars.last() != 'condition') {
        error("--design must end with 'condition', e.g. '~ batch + condition' (got '${design}')")
    }
    def missing = vars.findAll { v -> !rows[0].containsKey(v) }
    if (missing) {
        error("--design uses column(s) not in the samplesheet: ${missing.join(', ')}")
    }
    def empty = rows.findAll { r -> vars.any { v -> !r[v] } }.collect { r -> r.sample }
    if (empty) {
        error("Samplesheet has empty design values for sample(s): ${empty.join(', ')}")
    }
}

def checkContrasts(List levels) {
    if (levels.size() < 2) {
        error("The samplesheet needs at least two conditions (found: ${levels.join(', ')})")
    }
    if (!params.contrasts) {
        return
    }
    params.contrasts.toString().tokenize(',').collect { c -> c.trim() }.each { c ->
        def parts = c.split('_vs_')
        if (parts.size() != 2 || !(parts[0] in levels) || !(parts[1] in levels)) {
            error("Contrast '${c}' is not valid. Write LEVEL_vs_LEVEL using: ${levels.join(', ')}")
        }
    }
}

// ---------------------------------------------------------------------------
// Reference helpers
// ---------------------------------------------------------------------------

def genomeCache() {
    return params.genome_cache ?: "${System.getenv('HOME')}/.rnaseq-refs"
}

def ensemblUrls(String genome, release) {
    def species = [GRCh38: 'homo_sapiens', GRCm39: 'mus_musculus'][genome]
    def prefix  = [GRCh38: 'Homo_sapiens', GRCm39: 'Mus_musculus'][genome]
    def base    = "https://ftp.ensembl.org/pub/release-${release}"
    return [
        fasta: "${base}/fasta/${species}/dna/${prefix}.${genome}.dna.primary_assembly.fa.gz",
        gtf  : "${base}/gtf/${species}/${prefix}.${genome}.${release}.gtf.gz"
    ]
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {
    checkParams()
    def rows = readSamplesheet(params.input)
    checkDesign(rows)
    checkContrasts(rows.collect { r -> r.condition }.unique())

    // One item per sample: [ [id, single_end], [fastq_1, (fastq_2)] ]
    def sheet_dir = file(params.input).parent
    def reads = channel.fromList(rows).map { r ->
        def files = [r.fastq_1, r.fastq_2].findAll { p -> p }.collect { p -> resolvePath(p, sheet_dir) }
        [[id: r.sample, single_end: files.size() == 1], files]
    }

    // ---- Reference genome ------------------------------------------------
    // Ensembl downloads live in the shared cache; indexes for custom FASTA/GTF live in <outdir>/reference.
    def cache_dir = params.gtf
        ? "${params.outdir}/reference"
        : "${genomeCache()}/${params.genome}/ensembl_${params.ensembl_release}"
    def fasta = channel.empty()
    def gtf = channel.empty()
    if (params.gtf) {
        gtf = channel.value(file(params.gtf, checkIfExists: true))
        if (params.fasta) {
            fasta = channel.value(file(params.fasta, checkIfExists: true))
        }
    } else {
        def urls = ensemblUrls(params.genome, params.ensembl_release)
        DOWNLOAD_REFS(urls.fasta, urls.gtf, cache_dir)
        fasta = DOWNLOAD_REFS.out.fasta
        gtf = DOWNLOAD_REFS.out.gtf
    }
    def index = channel.empty()
    if (params.star_index) {
        index = channel.value(file(params.star_index, checkIfExists: true))
    } else {
        STAR_INDEX(fasta, gtf, cache_dir)
        index = STAR_INDEX.out.index
    }

    // ---- QC, trimming, alignment ----------------------------------------
    FASTQC(reads)
    def trimmed = reads
    def fastp_json = channel.empty()
    if (!params.skip_trimming) {
        FASTP(reads)
        trimmed = FASTP.out.reads
        fastp_json = FASTP.out.json
    }
    STAR_ALIGN(trimmed, index)

    // ---- Strandedness: auto-detect from STAR gene counts, or use --strandedness
    def strandedness = channel.value(params.strandedness)
    def strand_report = channel.value(file("${projectDir}/assets/NO_FILE"))
    if (params.strandedness == 'auto') {
        INFER_STRAND(STAR_ALIGN.out.gene_counts.collect())
        strandedness = INFER_STRAND.out.strand.map { f -> f.text.trim() }
        strand_report = INFER_STRAND.out.mqc
    }

    // ---- Counting and QC summaries --------------------------------------
    FEATURECOUNTS(STAR_ALIGN.out.bam, gtf, strandedness)
    MERGE_COUNTS(FEATURECOUNTS.out.counts.collect())
    GENE_NAMES(gtf)
    QC_SUMMARY(fastp_json.mix(STAR_ALIGN.out.log, FEATURECOUNTS.out.summary).collect())
    MULTIQC(
        FASTQC.out.zip
            .mix(fastp_json, STAR_ALIGN.out.log, FEATURECOUNTS.out.summary, strand_report)
            .collect()
    )

    // ---- Run records: samplesheet copy + tool versions (pipeline_info/) ---
    def samplesheet_copy = channel.fromPath(params.input)
        .collectFile(name: 'samplesheet.csv', storeDir: "${params.outdir}/pipeline_info")
    def version_lines = ["nextflow\t${workflow.nextflow.version}", "pipeline\t${workflow.manifest.version}"] +
        params.containers.collect { k, v -> "${k}\t${v}" }
    def versions = channel.fromList(version_lines)
        .collectFile(name: 'versions.tsv', newLine: true, sort: false, storeDir: "${params.outdir}/pipeline_info")
}
