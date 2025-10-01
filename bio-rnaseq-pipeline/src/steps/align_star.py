from pathlib import Path
from src.utils.shell import run_cmd

def run_star(samples_df, cfg):
    trimmed = Path(cfg["outdir"]).joinpath("trimmed")
    star_out = Path(cfg["outdir"]).joinpath("star")
    idx = cfg["reference"]["star_index"]
    sparams = cfg["star"]
    for _, s in samples_df.iterrows():
        prefix = star_out / s["sample"]
        if s["mode"] == "paired":
            r1 = trimmed/f"{s['sample']}_R1.trim.fastq.gz"
            r2 = trimmed/f"{s['sample']}_R2.trim.fastq.gz"
            reads = f"--readFilesIn {r1} {r2} --readFilesCommand zcat"
        else:
            r1 = trimmed/f"{s['sample']}_R1.trim.fastq.gz"
            reads = f"--readFilesIn {r1} --readFilesCommand zcat"
        cmd = (
            f"STAR --runThreadN {cfg['threads']} --genomeDir {idx} {reads} "
            f"--outFileNamePrefix {prefix}. "
            f"--outFilterMismatchNmax {sparams['out_filter_mismatch_nmax']} "
            f"--outSAMtype {sparams['out_sam_type']}"
        )
        run_cmd(cmd)
