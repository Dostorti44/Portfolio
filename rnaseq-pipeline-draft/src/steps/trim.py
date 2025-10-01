from pathlib import Path
from src.utils.shell import run_cmd

def run_fastp(samples_df, cfg):
    trimmed_dir = Path(cfg["outdir"]).joinpath("trimmed")
    p = cfg["fastp"]
    for _, s in samples_df.iterrows():
        if s["mode"] == "paired":
            o1 = trimmed_dir/f"{s['sample']}_R1.trim.fastq.gz"
            o2 = trimmed_dir/f"{s['sample']}_R2.trim.fastq.gz"
            html = trimmed_dir/f"{s['sample']}.fastp.html"
            json = trimmed_dir/f"{s['sample']}.fastp.json"
            cmd = (
                f"fastp -i {s['r1']} -I {s['r2']} -o {o1} -O {o2} "
                f"--trim_front1 {p['trim_front1']} --trim_tail1 {p['trim_tail1']} "
                f"--cut_mean_quality {p['cut_mean_quality']} --length_required {p['length_required']} "
                f"--thread {cfg['threads']} --html {html} --json {json}"
            )
        else:
            o1 = trimmed_dir/f"{s['sample']}_R1.trim.fastq.gz"
            html = trimmed_dir/f"{s['sample']}.fastp.html"
            json = trimmed_dir/f"{s['sample']}.fastp.json"
            cmd = (
                f"fastp -i {s['r1']} -o {o1} "
                f"--trim_front1 {p['trim_front1']} --trim_tail1 {p['trim_tail1']} "
                f"--cut_mean_quality {p['cut_mean_quality']} --length_required {p['length_required']} "
                f"--thread {cfg['threads']} --html {html} --json {json}"
            )
        run_cmd(cmd)
