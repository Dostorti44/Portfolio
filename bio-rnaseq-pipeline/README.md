# RNA‑seq Pipeline (STAR → featureCounts → DESeq2)

Reproducible RNA‑seq analysis from raw FASTQs to differential expression.
Python-first CLI with R/DESeq2 via rpy2 for statistics.

## Quick start
```bash
mamba env create -f environment.yml
conda activate rnaseq-pipeline
python -m src.main all
```

## Inputs
- `samples/samples.csv`: sample manifest (`sample,condition,r1,r2` for paired; `sample,condition,r1` for single).
- `refs/`: `genome.fa` and `genes.gtf` (do **not** commit these).
- `config/config.yaml`: tool parameters and DE design.

## Outputs
- `results/fastqc`, `results/reports/multiqc_report.html`
- `results/trimmed`: trimmed FASTQs
- `results/star`: sorted BAMs + logs
- `results/counts/gene_counts.txt`: gene x sample matrix
- `results/de`: DESeq2 results + PCA
- `results/reports/volcano.png`, `top50_DE_genes.csv`

## Rationale
- **STAR** for spliced alignment; **featureCounts** for gene counts; **DESeq2** for differential expression.

## Reproducibility
- Conda environment pinned in `environment.yml`.
- Config-driven parameters in `config/config.yaml`.

## Citation
Please cite STAR, featureCounts, FastQC, MultiQC, fastp, and DESeq2 when publishing results.
