from pathlib import Path
from src.utils.shell import run_cmd
import pandas as pd

def run_featurecounts(samples_df, cfg):
    star_out = Path(cfg["outdir"]).joinpath("star")
    counts_dir = Path(cfg["outdir"]).joinpath("counts")
    gtf = cfg["reference"]["gtf"]
    paired = cfg["featurecounts"]["is_paired_end"]
    strand = cfg["featurecounts"]["strand"]

    bam_list = [str(star_out/f"{s['sample']}.Aligned.sortedByCoord.out.bam") for _, s in samples_df.iterrows()]
    bam_str = " ".join(bam_list)
    out = counts_dir/"gene_counts.txt"

    paired_flag = "-p" if paired else ""
    cmd = (
        f"featureCounts -T {cfg['threads']} -a {gtf} -o {out} {paired_flag} -s {strand} "
        f"-g {cfg['featurecounts']['gtf_attr']} {bam_str}"
    )
    run_cmd(cmd)

    samples_df[['sample','condition']].to_csv(counts_dir/'sample_table.csv', index=False)
