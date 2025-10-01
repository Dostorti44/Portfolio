from pathlib import Path
from src.utils.shell import run_cmd

def build_star_index(cfg):
    idx = Path(cfg["reference"]["star_index"])
    if list(idx.glob("*")):
        return  # already built
    fa = cfg["reference"]["fasta"]
    gtf = cfg["reference"]["gtf"]
    sj = cfg["star"]["sjdb_overhang"]
    threads = cfg["threads"]
    cmd = (
        f"STAR --runThreadN {threads} --runMode genomeGenerate "
        f"--genomeDir {idx} --genomeFastaFiles {fa} --sjdbGTFfile {gtf} "
        f"--sjdbOverhang {sj}"
    )
    run_cmd(cmd)
