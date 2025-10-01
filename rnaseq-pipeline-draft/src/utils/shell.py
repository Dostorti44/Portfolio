import subprocess as sp
from rich import print

def run_cmd(cmd: str, cwd=None):
    print(f"[dim]$ {cmd}[/]")
    res = sp.run(cmd, shell=True, cwd=cwd)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}")
