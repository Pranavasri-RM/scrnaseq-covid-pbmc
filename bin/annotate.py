#!/usr/bin/env python3
"""Remove artifact cells, name clusters from a labels file and split the NK / CD8 T cluster with marker scores."""

import argparse
import os
import warnings

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)   # celltypist uses a deprecated scanpy attribute on import
    import celltypist
    from celltypist import models

# Marker sets for splitting the NK / CD8 T cluster (several genes each, since single genes are often not detected)
T_GENES  = ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B"]
NK_GENES = ["KLRF1", "NCAM1", "KLRD1", "TYROBP", "FCER1G"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="clustered .h5ad")
    p.add_argument("--labels", required=True,
                   help="CSV with columns cluster, cell_type, marker (REMOVE drops a cluster)")
    p.add_argument("--model", default="Immune_All_Low.pkl", help="CellTypist model name or path to a .pkl file")
    p.add_argument("--out", required=True, help="annotated .h5ad")
    p.add_argument("--counts-out", required=True, help="CSV of cells per sample and cell type")
    return p.parse_args()


def gene_ids(ad, symbols):
    """Ensembl IDs for gene names (first match), skipping names that aren't in the data."""
    first = ad.var.reset_index().drop_duplicates("gene_name").set_index("gene_name").iloc[:, 0]
    return [first[s] for s in symbols if s in first.index]


def check_markers(ad, labels_df):
    """Each labelled cluster must have the highest mean expression of its marker gene."""
    clusters = ad.obs["leiden_0.5"].astype(str)
    for cluster, row in labels_df.dropna(subset=["marker"]).iterrows():
        gene = ad.var_names[ad.var["gene_name"] == row["marker"]]
        assert len(gene) == 1, f"marker {row['marker']} not found once in the genes"
        expr = pd.Series(ad[:, gene[0]].X.toarray().ravel(), index=ad.obs_names)
        top = expr.groupby(clusters).mean().idxmax()
        assert top == cluster, (f"{row['marker']} ({row['cell_type']}) peaks in cluster {top}, "
                                f"not cluster {cluster}: cluster labels don't match this clustering")
    print("Marker check passed: every marker peaks in its labelled cluster")


def split_nk_cd8(ad, mixed_label="NK / CD8 T"):
    """Label cells of the mixed cluster as CD8 T or NK by T vs NK marker scores, smoothed over neighbors."""
    mixed = (ad.obs["cell_type"] == mixed_label).values
    sub = ad[mixed].copy()
    sc.tl.score_genes(sub, gene_ids(ad, T_GENES), score_name="t_score", random_state=0)
    sc.tl.score_genes(sub, gene_ids(ad, NK_GENES), score_name="nk_score", random_state=0)

    # Average each cell's scores with its neighbors (graph from cluster.py, restricted to this cluster)
    W = ad.obsp["connectivities"][mixed][:, mixed] + sp.identity(mixed.sum())
    W = sp.diags(1 / np.asarray(W.sum(axis=1)).ravel()) @ W
    t_s  = W @ sub.obs["t_score"].values
    nk_s = W @ sub.obs["nk_score"].values

    ad.obs["t_score"], ad.obs["nk_score"] = np.nan, np.nan
    ad.obs.loc[mixed, "t_score"], ad.obs.loc[mixed, "nk_score"] = t_s, nk_s
    ad.obs["cell_type"] = ad.obs["cell_type"].astype(str)
    ad.obs.loc[mixed, "cell_type"] = np.where(t_s > nk_s, "CD8 T", "NK")


def main():
    args = parse_args()
    anndata.settings.allow_write_nullable_strings = True

    ad = sc.read_h5ad(args.input)
    labels_df = pd.read_csv(args.labels, dtype=str).set_index("cluster")
    labels = labels_df["cell_type"]
    clusters = ad.obs["leiden_0.5"].astype(str)
    assert set(clusters) <= set(labels.index), f"clusters without a label: {sorted(set(clusters) - set(labels.index))}"
    check_markers(ad, labels_df)

    # Artifact cells: nuclear-RNA-like (high genes per UMI and very low MT), plus clusters labelled REMOVE
    genes_per_umi = ad.obs["n_genes_by_counts"] / ad.obs["total_counts"]
    artifact = (genes_per_umi > 0.75) & (ad.obs["pct_counts_mt"] < 2)
    remove = artifact | (clusters.map(labels) == "REMOVE")
    ad = ad[~remove].copy()
    print(f"Removed {remove.sum()} artifact cells, {ad.n_obs} left")

    ad.obs["cell_type"] = ad.obs["leiden_0.5"].astype(str).map(labels)

    # CellTypist labels, kept as an independent comparison (needs gene names instead of Ensembl IDs)
    if not os.path.exists(args.model):
        models.download_models(model=[args.model])
    ct_input = anndata.AnnData(X=ad.X, obs=ad.obs[[]], var=pd.DataFrame(index=ad.var["gene_name"].values))
    ct_input.var_names_make_unique()
    pred = celltypist.annotate(ct_input, model=args.model, majority_voting=True,
                               over_clustering=ad.obs["leiden_1.0"].values)
    ad.obs["celltypist"] = pred.predicted_labels["majority_voting"].values
    ad.obs["celltypist_cell"] = pred.predicted_labels["predicted_labels"].values

    split_nk_cd8(ad)
    ad.obs["cell_type"] = ad.obs["cell_type"].astype("category")

    ad.write_h5ad(args.out)
    pd.crosstab(ad.obs["sample"], ad.obs["cell_type"]).to_csv(args.counts_out)
    print(ad.obs["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()