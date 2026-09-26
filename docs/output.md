# Outputs

```
results/
├── fastqc/               FastQC reports for raw reads (<sample>_fastqc.html)
├── fastp/                Trimming reports (<sample>.fastp.html / .json)
├── star/                 <sample>.Log.final.out (mapping stats), <sample>.ReadsPerGene.out.tab,
│                         strandedness.txt + strandedness_mqc.tsv (auto mode), BAMs with --save_bam
│                         (with --bam_dir the BAMs go to that folder instead)
├── featurecounts/        counts.tsv (genes x samples, raw counts), gene_names.tsv,
│   └── per_sample/       featureCounts output and summary per sample
├── qc/qc_summary.tsv     Reads, % retained, % uniquely mapped, % assigned per sample
├── multiqc/              multiqc_report.html: all upstream QC in one page
├── report/               report.html + tables/ (CSV per contrast); not created with --skip_downstream
├── reference/            STAR index, only when --fasta/--gtf were given
└── pipeline_info/        samplesheet.csv, versions.tsv, timeline/trace/execution reports
```

`counts.tsv` holds **raw** counts: use it for DESeq2/edgeR, not for plotting expression directly.

## Result tables (`report/tables/`)

`<A>_vs_<B>_deseq2.csv`: `gene_id, gene_name, baseMean, log2FC, log2FC_shrunk, pvalue, padj, significant`
`<A>_vs_<B>_edger.csv`: `gene_id, gene_name, logCPM, log2FC, pvalue, padj, significant`
`<A>_vs_<B>_go_up.csv` / `_go_down.csv`: enriched GO terms (only written when terms are found)

Positive log2FC = higher in A. `significant` = padj < `padj_cutoff` and |log2FC| > `lfc_cutoff`
(unshrunk log2FC for both methods). `log2FC_shrunk` (apeglm) is better for ranking genes and
for plots because noisy low-count genes are pulled towards 0.

## Reading the report

**Warnings and checks.** Read first. Lists anything that changed the analysis: a design that couldn't be
estimated (replaced by `~ condition`), groups with < 2 replicates, low gene assignment.

**1. Run summary.** Genome, Ensembl release, strandedness used, design, contrasts, cutoffs, tool versions.

**2. Upstream QC.** Per-sample table and bars. Uniquely mapped is typically 70–90% for good human/mouse
data; below 60% is flagged. Assigned to genes below 50% is flagged, usually a strandedness or annotation issue.

**3. Count QC.** Library sizes should be of the same order across samples. The log-CPM boxes
should look similar; one sample far off is worth checking.

**4. Sample relationships.**
- *PCA*: each dot is a sample. Replicates of a condition should sit together; if samples split by batch
  instead, batch is a strong effect and belongs in `--design`.
- *Scree*: how much variation each PC captures.
- *Distance heatmap*: darker = more similar. Look for outliers.
- *PCA after removing covariates*: shows what the data look like once batch is accounted for. For display only.

**5. Differential expression** (one tab per contrast).
- *MA plot*: expression (x) vs fold change (y). Red points are significant.
- *Volcano*: fold change (x) vs −log10 p-value (y). Top genes labelled.
- *Heatmap*: top significant genes, scaled per gene (z-score); columns should group by condition.
- *Tables*: searchable, top 1000 genes; full tables in `tables/`.

**6. DESeq2 vs edgeR.** Genes significant in both methods are the most robust. Fold changes
normally correlate strongly (Spearman > 0.9).

**7. GO enrichment.** Over-represented Biological Process / Molecular Function / Cellular Component
terms among significant DESeq2 genes, up and down separately. Dot size = number of genes,
colour = adjusted p-value. Fewer than 5 significant genes: skipped.
