from pathlib import Path
import pandas as pd
import yaml

REQUIRED_COLUMNS_PAIRED = ["sample", "condition", "r1", "r2"]
REQUIRED_COLUMNS_SINGLE = ["sample", "condition", "r1"]

def read_config(path: str) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["outdir"] = str(Path(cfg["outdir"]).resolve())
    return cfg

def read_samples(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    cols = set(df.columns)
    if set(REQUIRED_COLUMNS_PAIRED).issubset(cols):
        mode = "paired"
    elif set(REQUIRED_COLUMNS_SINGLE).issubset(cols):
        mode = "single"
    else:
        raise ValueError("samples.csv missing required columns")
    df["mode"] = mode
    for col in ["r1", "r2"]:
        if col in df.columns:
            missing = [p for p in df[col].dropna() if not Path(p).exists()]
            if missing:
                raise FileNotFoundError(f"Missing FASTQs: {missing[:3]} ...")
    return df

def ensure_dirs(cfg: dict):
    for sub in ["fastqc", "trimmed", "star", "counts", "de", "reports"]:
        Path(cfg["outdir"]).joinpath(sub).mkdir(parents=True, exist_ok=True)
    Path(cfg["reference"]["star_index"]).mkdir(parents=True, exist_ok=True)
