from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def build_report(cfg: dict):
    de_dir = Path(cfg["outdir"]).joinpath("de")
    reports = Path(cfg["outdir"]).joinpath("reports")
    reports.mkdir(parents=True, exist_ok=True)

    res = pd.read_csv(de_dir/"deseq2_results.csv", index_col=0)
    padj_thr = cfg["report"]["volcano_padj"]
    lfc_thr = cfg["report"]["volcano_lfc"]

    x = res["log2FoldChange"]
    y = -np.log10(res["padj"].fillna(1.0))

    sig = (res["padj"] < padj_thr) & (res["log2FoldChange"].abs() >= lfc_thr)

    plt.figure()
    plt.scatter(x[~sig], y[~sig], s=6, alpha=0.5)
    plt.scatter(x[sig], y[sig], s=6, alpha=0.8)
    plt.xlabel("log2 fold-change")
    plt.ylabel("-log10 adjusted p-value")
    plt.title("Volcano plot")
    plt.savefig(reports/"volcano.png", dpi=200, bbox_inches="tight")
    plt.close()

    top = res.dropna(subset=["padj"]).sort_values("padj").head(50)
    top.to_csv(reports/"top50_DE_genes.csv")
