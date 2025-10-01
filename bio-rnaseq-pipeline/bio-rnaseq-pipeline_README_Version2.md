# RNA-seq Analysis Pipeline

A reproducible, end-to-end RNA-seq pipeline from raw FASTQ files to differential gene expression results. This project provides a command-line interface (CLI) in Python, integrating R/DESeq2 for statistical analysis. The workflow is designed to be simple to use, customizable, and robust for both small and large-scale RNA-seq experiments.

## Overview

This pipeline automates standard RNA-seq analysis steps:

1. **Quality control** (FastQC, MultiQC)
2. **Read trimming** (fastp)
3. **Genome indexing** (STAR)
4. **Read alignment** (STAR)
5. **Gene counting** (featureCounts)
6. **Differential expression** (DESeq2 via rpy2)
7. **Reporting** (PCA, volcano plot, top genes)

It is built for reproducibility, with all parameters configurable and software managed via Conda.

## Installation

Clone this repository and set up the environment:

```bash
git clone https://github.com/Dostorti44/Portfolio.git
cd Portfolio/bio-rnaseq-pipeline
mamba env create -f environment.yml
conda activate rnaseq-pipeline
```

## Usage

Run the entire workflow with the default configuration and sample files:

```bash
python -m src.main all
```

Or specify custom paths:

```bash
python -m src.main all --config-path config/config.yaml --samples-csv samples/samples.csv
```

## Input Files

| File                           | Description                                                                                           |
|---------------------------------|-------------------------------------------------------------------------------------------------------|
| `samples/samples.csv`           | Sample manifest. Columns: `sample`, `condition`, `r1`, `r2` (paired) or `sample`, `condition`, `r1` (single) |
| `refs/genome.fa`                | Reference genome (FASTA format)                                                                       |
| `refs/genes.gtf`                | Gene annotation (GTF format)                                                                          |
| `config/config.yaml`            | Pipeline configuration: tool parameters, DE design, output paths                                      |
| `config/contrasts.csv`          | Contrasts for differential analysis (e.g. Treated vs Control)                                         |

> **Note:** Do not commit large reference files (`genome.fa`, `genes.gtf`) to source control.

## Output Files

| Directory / File                       | Contents                                                               |
|-----------------------------------------|------------------------------------------------------------------------|
| `results/fastqc`                       | FastQC reports for raw reads                                           |
| `results/reports/multiqc_report.html`   | MultiQC summary                                                        |
| `results/trimmed`                      | Trimmed FASTQ files                                                    |
| `results/star`                         | Aligned, sorted BAM files and STAR logs                                |
| `results/counts/gene_counts.txt`        | Gene count matrix (genes × samples)                                    |
| `results/de/deseq2_results.csv`         | DESeq2 results (all genes, statistics)                                 |
| `results/de/PCA.pdf`                    | Principal Component Analysis plot                                      |
| `results/reports/volcano.png`           | Volcano plot (DE results visualization)                                |
| `results/reports/top50_DE_genes.csv`    | Top 50 differentially expressed genes                                  |

## Configuration

All major parameters are controlled via `config/config.yaml`, including:

- Thread count
- Reference file paths
- Trimming and alignment settings
- featureCounts options
- DESeq2 design formula and reference condition
- Reporting thresholds

Example config fragment:

```yaml
project: "rnaseq-demo"
threads: 8
outdir: "results"
reference:
  fasta: "refs/genome.fa"
  gtf: "refs/genes.gtf"
  star_index: "refs/star_index"
# ... more options ...
```

## Pipeline Workflow

1. **Quality Control**: Runs FastQC on input FASTQ files; MultiQC consolidates reports.
2. **Trimming**: Uses fastp to trim adapters and low-quality bases.
3. **Genome Indexing**: STAR index built from `genome.fa` and `genes.gtf` (skipped if present).
4. **Alignment**: STAR aligns trimmed reads, outputs sorted BAM files.
5. **Counting**: featureCounts quantifies gene-level read counts.
6. **Differential Expression**: DESeq2 (via rpy2) analyzes counts for differential expression.
7. **Reporting**: Outputs volcano plot, PCA, and top DE genes.

## Extending & Customizing

- Adjust parameters in `config/config.yaml` for specific experiment needs.
- Add new analysis steps or modify scripts in `src/steps/`.
- The modular structure makes it easy to swap tools or add custom logic.

## Troubleshooting

- Ensure all FASTQ, reference, and annotation files are present and paths are correct.
- Check Conda environment for required packages.
- Output and logs are saved in the `results/` directory for review.

## Citation & Credits

If you use this pipeline, please cite the following tools in your publication:

- STAR
- featureCounts
- FastQC
- MultiQC
- fastp
- DESeq2

## License

This pipeline is released under the MIT License. See [LICENSE](LICENSE) for details.

---

For questions, suggestions, or contributions, please open an issue or pull request.