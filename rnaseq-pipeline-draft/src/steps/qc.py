from pathlib import Path
from src.utils.shell import run_cmd

def run_fastqc_multiqc(samples_df, cfg):
    out_fastqc = Path(cfg["outdir"]).joinpath("fastqc")
    fastqs = []
    for _, row in samples_df.iterrows():
        fastqs.append(row["r1"])
        if row["mode"] == "paired":
            fastqs.append(row["r2"])
    fastq_str = " ".join(fastqs)
    run_cmd(f"fastqc -t {cfg['threads']} -o {out_fastqc} {fastq_str}")
    run_cmd(f"multiqc {Path(cfg['outdir'])} -o {Path(cfg['outdir']).joinpath('reports')}")
