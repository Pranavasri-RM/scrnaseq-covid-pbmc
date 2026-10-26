#!/usr/bin/env python3
"""Combine per-sample QC outputs, drop rarely detected genes and keep raw counts in a layer."""

import argparse
from pathlib import Path

import anndata
import pandas as pd
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", nargs="+", required=True, help="per-sample QC .h5ad files named <sample>.h5ad")
    p.add_argument("--samplesheet", required=True, help="sample sheet; samples are combined in its order")
    p.add_argument("--min-cells", type=int, default=10, help="keep genes detected in at least this many cells")
    p.add_argument("--out", required=True, help="combined .h5ad")
    return p.parse_args()


def main():
    args = parse_args()
    anndata.settings.allow_write_nullable_strings = True

    # Combine in sample sheet order, so the result is the same on every run
    order = pd.read_csv(args.samplesheet)["sample"].tolist()
    files = {Path(f).stem: f for f in args.inputs}
    assert set(files) == set(order), f"inputs don't match the sample sheet: {sorted(set(files) ^ set(order))}"
    adatas = [sc.read_h5ad(files[s]) for s in order]

    ad = sc.concat(adatas, merge="same")   # keep var columns that are identical across samples

    n_genes = ad.n_vars
    sc.pp.filter_genes(ad, min_cells=args.min_cells)
    ad.layers["counts"] = ad.X.copy()      # raw counts kept before any normalization

    ad.write_h5ad(args.out)
    print(f"Combined: {ad.n_obs} cells from {len(adatas)} samples")
    print(f"Genes: {ad.n_vars} of {n_genes} detected in at least {args.min_cells} cells")


if __name__ == "__main__":
    main()