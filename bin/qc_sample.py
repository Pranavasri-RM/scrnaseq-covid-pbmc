#!/usr/bin/env python3
"""Per-sample QC: add gene names and QC metrics, apply the saved cutoffs, score doublets."""

import argparse

import anndata
import pandas as pd
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts-dir", required=True, help="kb count output folder (counts_unfiltered/)")
    p.add_argument("--sample", required=True, help="sample name, as in the sample sheet")
    p.add_argument("--donor", required=True)
    p.add_argument("--condition", required=True)
    p.add_argument("--cutoffs", required=True, help="CSV with columns sample, min_umi, max_umi, max_pct_mt")
    p.add_argument("--out", required=True, help="output .h5ad")
    return p.parse_args()


def load_counts(counts_dir):
    """Read the kb count matrix and add gene names, checking they line up with the matrix columns."""
    ad = sc.read_h5ad(f"{counts_dir}/adata.h5ad")
    ids   = pd.read_csv(f"{counts_dir}/cells_x_genes.genes.txt",       header=None, skip_blank_lines=False)[0]
    names = pd.read_csv(f"{counts_dir}/cells_x_genes.genes.names.txt", header=None, skip_blank_lines=False)[0]
    assert len(ids) == len(names) == ad.n_vars, "gene count mismatch"
    assert (ids.values == ad.var_names).all(), "gene order differs from the matrix columns"
    ad.var["gene_name"] = names.fillna(ids).values   # unnamed genes keep their ID
    return ad


def main():
    args = parse_args()
    anndata.settings.allow_write_nullable_strings = True

    ad = load_counts(args.counts_dir)

    # QC metrics
    ad.var["mt"] = ad.var["gene_name"].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Sample information on every cell; barcodes made unique across samples
    ad.obs["sample"], ad.obs["donor"], ad.obs["condition"] = args.sample, args.donor, args.condition
    ad.obs_names = args.sample + "_" + ad.obs_names

    # Cutoffs chosen in notebooks/01_qc_exploration.ipynb
    cut = pd.read_csv(args.cutoffs, index_col="sample").loc[args.sample]
    keep = ((ad.obs["total_counts"] >= cut["min_umi"]) &
            (ad.obs["total_counts"] <= cut["max_umi"]) &
            (ad.obs["pct_counts_mt"] <= cut["max_pct_mt"]))
    n_before = ad.n_obs
    ad = ad[keep].copy()

    # Doublet scores within this sample (cutoff not applied, see the QC notebook)
    sc.pp.scrublet(ad, expected_doublet_rate=0.05, random_state=0)

    ad.write_h5ad(args.out)
    print(f"{args.sample}: kept {ad.n_obs} of {n_before} barcodes")


if __name__ == "__main__":
    main()