from pathlib import Path
import rpy2.robjects as ro

R = ro.r

R(
    """

    suppressMessages({
      library(DESeq2)
      library(tidyverse)
    })

    """

)

def run_deseq2(cfg: dict):
    counts_fp = Path(cfg["outdir"]).joinpath("counts", "gene_counts.txt")
    sample_table_fp = Path(cfg["outdir"]).joinpath("counts", "sample_table.csv")
    de_dir = Path(cfg["outdir"]).joinpath("de")
    de_dir.mkdir(parents=True, exist_ok=True)

    R.assign("counts_fp", str(counts_fp))
    R.assign("sample_table_fp", str(sample_table_fp))
    R.assign("design_formula", cfg["design"]["formula"])
    R.assign("ref_level", cfg["design"]["ref_level"])
    R.assign("de_dir", str(de_dir))

    R(
        """

        fc <- read.delim(counts_fp, comment.char="#")
        rownames(fc) <- fc$Geneid
        counts <- as.matrix(fc[,7:ncol(fc)])

        coldata <- read.csv(sample_table_fp)
        stopifnot(all(colnames(counts) == coldata$sample))
        coldata$condition <- factor(coldata$condition)
        coldata$condition <- relevel(coldata$condition, ref=ref_level)

        dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=as.formula(design_formula))
        dds <- dds[rowSums(counts(dds)) > 10, ]
        dds <- DESeq(dds)
        res <- results(dds)
        resOrdered <- res[order(res$padj),]

        write.csv(as.data.frame(resOrdered), file=file.path(de_dir, "deseq2_results.csv"))
        vsd <- vst(dds)
        write.csv(as.data.frame(assay(vsd)), file=file.path(de_dir, "vsd_matrix.csv"))

        pdf(file.path(de_dir, "PCA.pdf"))
        plotPCA(vsd, intgroup="condition")
        dev.off()
        """

    )
