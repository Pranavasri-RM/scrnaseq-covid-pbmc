#!/usr/bin/env python3
"""Remove artifact cells, name clusters from a labels file and split NK / CD8 T cells with CellTypist."""

import argparse
import os
import warnings

import anndata
import pandas as pd
import scanpy as sc

with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)   # celltypist uses a deprecated scanpy attribute on import
    import celltypist
    from celltypist import models


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="clustered .h5ad")
    p.add_argument("--labels", required=True, help="CSV with columns cluster, cell_type (REMOVE drops a cluster)")
    p.add_argument("--model", default="Immune_All_Low.pkl", help="CellTypist model name or path to a .pkl file")
    p.add_argument("--out", required=True, help="annotated .h5ad")
    p.add_argument("--counts-out", required=True, help="CSV of cells per sample and cell type")
    return p.parse_args()


def main():
    args = parse_args()
    anndata.settings.allow_write_nullable_strings = True

    ad = sc.read_h5ad(args.input)
    labels = pd.read_csv(args.labels, dtype=str).set_index("cluster")["cell_type"]
    clusters = ad.obs["leiden_0.5"].astype(str)
    assert set(clusters) <= set(labels.index), f"clusters without a label: {sorted(set(clusters) - set(labels.index))}"

    # Artifact cells: nuclear-RNA-like (high genes per UMI and very low MT), plus clusters labelled REMOVE
    genes_per_umi = ad.obs["n_genes_by_counts"] / ad.obs["total_counts"]
    artifact = (genes_per_umi > 0.75) & (ad.obs["pct_counts_mt"] < 2)
    remove = artifact | (clusters.map(labels) == "REMOVE")
    ad = ad[~remove].copy()
    print(f"Removed {remove.sum()} artifact cells, {ad.n_obs} left")

    ad.obs["cell_type"] = ad.obs["leiden_0.5"].astype(str).map(labels)

    # CellTypist (majority voting over leiden_1.0); needs gene names instead of Ensembl IDs
    if not os.path.exists(args.model):
        models.download_models(model=[args.model])
    ct_input = anndata.AnnData(X=ad.X, obs=ad.obs[[]], var=pd.DataFrame(index=ad.var["gene_name"].values))
    ct_input.var_names_make_unique()
    pred = celltypist.annotate(ct_input, model=args.model, majority_voting=True,
                               over_clustering=ad.obs["leiden_1.0"].values)
    ad.obs["celltypist"] = pred.predicted_labels["majority_voting"].values
    ad.obs["celltypist_cell"] = pred.predicted_labels["predicted_labels"].values

    # Split the NK / CD8 T cluster with CellTypist's majority-vote labels
    ct = ad.obs["celltypist"].astype(str)
    mixed = ad.obs["cell_type"] == "NK / CD8 T"
    ad.obs.loc[mixed & ct.str.contains("NK"), "cell_type"] = "NK"
    ad.obs.loc[mixed & ct.str.contains("helper"), "cell_type"] = "CD4 T"
    ad.obs.loc[mixed & ~ct.str.contains("NK|helper"), "cell_type"] = "CD8 T"
    ad.obs["cell_type"] = ad.obs["cell_type"].astype("category")

    ad.write_h5ad(args.out)
    pd.crosstab(ad.obs["sample"], ad.obs["cell_type"]).to_csv(args.counts_out)
    print(ad.obs["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()