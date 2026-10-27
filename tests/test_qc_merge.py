import subprocess
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "qc_merge.py"


def write_sample(path, name, n_cells):
    ad = anndata.AnnData(X=np.ones((n_cells, 3)))
    ad.obs_names = [f"{name}_cell{i}" for i in range(n_cells)]
    ad.obs["sample"] = name
    ad.var_names = ["G1", "G2", "G3"]
    ad.var["gene_name"] = ["A", "B", "C"]
    ad.write_h5ad(path / f"{name}.h5ad")


def run_merge(tmp_path, inputs, order):
    pd.DataFrame({"sample": order}).to_csv(tmp_path / "sheet.csv", index=False)
    return subprocess.run([sys.executable, str(SCRIPT), "--inputs", *map(str, inputs),
                           "--samplesheet", str(tmp_path / "sheet.csv"), "--min-cells", "1",
                           "--out", str(tmp_path / "out.h5ad")], capture_output=True, text=True)


def test_merge_follows_sample_sheet_order(tmp_path):
    write_sample(tmp_path, "zeta", 2)
    write_sample(tmp_path, "alpha", 3)
    result = run_merge(tmp_path, [tmp_path / "alpha.h5ad", tmp_path / "zeta.h5ad"], order=["zeta", "alpha"])
    assert result.returncode == 0, result.stderr
    out = anndata.read_h5ad(tmp_path / "out.h5ad")
    assert list(out.obs["sample"]) == ["zeta"] * 2 + ["alpha"] * 3   # sheet order, not file-name order
    assert "counts" in out.layers


def test_merge_rejects_files_not_in_sample_sheet(tmp_path):
    # The bug from local testing: a stray combined file picked up by a glob
    write_sample(tmp_path, "alpha", 3)
    write_sample(tmp_path, "combined_qc", 3)
    result = run_merge(tmp_path, [tmp_path / "alpha.h5ad", tmp_path / "combined_qc.h5ad"], order=["alpha"])
    assert result.returncode != 0
    assert "don't match the sample sheet" in result.stderr