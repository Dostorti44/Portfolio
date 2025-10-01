# src/main.py
from pathlib import Path
import typer
from rich import print
from src.utils.io import read_config, read_samples, ensure_dirs
from src.steps.qc import run_fastqc_multiqc
from src.steps.trim import run_fastp
from src.steps.index_star import build_star_index
from src.steps.align_star import run_star
from src.steps.count_featurecounts import run_featurecounts
from src.steps.de_deseq2 import run_deseq2
from src.steps.report import build_report

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def all(config_path: str = "config/config.yaml", samples_csv: str = "samples/samples.csv"):
    cfg = read_config(config_path)
    samples = read_samples(samples_csv)
    ensure_dirs(cfg)

    print("[bold cyan]Step 1/6: QC (FastQC + MultiQC)...[/]")
    run_fastqc_multiqc(samples, cfg)

    print("[bold cyan]Step 2/6: Trimming (fastp)...[/]")
    run_fastp(samples, cfg)

    print("[bold cyan]Step 3/6: STAR index...[/]")
    build_star_index(cfg)

    print("[bold cyan]Step 4/6: STAR alignment...[/]")
    run_star(samples, cfg)

    print("[bold cyan]Step 5/6: featureCounts...[/]")
    run_featurecounts(samples, cfg)

    print("[bold cyan]Step 6/6: DESeq2 (rpy2) + report...[/]")
    run_deseq2(cfg)
    build_report(cfg)

    print("[bold green]Done. See results/ for outputs.[/]")

if __name__ == "__main__":
    app()
