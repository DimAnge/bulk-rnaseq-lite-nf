# R environment for the bulk-rnaseq-lite-nf HTML report.
# Build once:  docker build -t bulk-rnaseq-lite-report:1.0 docker/
# Bioconductor 3.22 pins package versions; the base image already has R, pandoc and rmarkdown.
FROM bioconductor/bioconductor_docker:RELEASE_3_22

RUN R -q -e 'pkgs <- c("DESeq2", "edgeR", "limma", "apeglm", "clusterProfiler", "enrichplot", \
                       "org.Hs.eg.db", "org.Mm.eg.db", "rmarkdown", "knitr", "ggplot2", "pheatmap", "DT"); \
             BiocManager::install(pkgs, ask = FALSE, update = FALSE); \
             missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
             if (length(missing) > 0) stop("Failed to install: ", paste(missing, collapse = ", "))'
