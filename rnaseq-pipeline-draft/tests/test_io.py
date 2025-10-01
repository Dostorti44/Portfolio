from src.utils.io import read_config
import tempfile, pathlib

def test_read_config():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "config.yaml"
        p.write_text("outdir: results\nreference: {fasta: refs/genome.fa, gtf: refs/genes.gtf, star_index: refs/star_index}")
        cfg = read_config(str(p))
        assert "outdir" in cfg
