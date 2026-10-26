#!/usr/bin/env python3
"""Normalize, select HVGs, run PCA, integrate samples with Harmony, then UMAP and Leiden clustering."""

import argparse

import anndata
import harmonypy
import numpy as np
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="combined QC .h5ad (raw counts in layers['counts'])")
    p.add_argument("--out", required=True, help="clustered .h5ad")
    p.add_argument("--n-hvgs", type=int, default=2000)
    p.add_argument("--n-pcs", type=int, default=30, help="PCs used for neighbors (chosen in notebook 02)")
    p.add_argument("--resolutions", type=float, nargs="+", default=[0.5, 1.0])
    return p.parse_args()


def run_harmony(ad, key):
    """Harmony on the PCA, called directly (scanpy's wrapper fails with harmonypy 2.x)."""
    ho = harmonypy.run_harmony(ad.obsm["X_pca"], ad.obs, key, random_state=0)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != ad.n_obs:      # older harmonypy returns PCs x cells
        Z = Z.T
    assert Z.shape == ad.obsm["X_pca"].shape
    return Z


def main():
    args = parse_args()
    anndata.settings.allow_write_nullable_strings = True

    ad = sc.read_h5ad(args.input)

    # Normalize: 10,000 UMIs per cell, then log(x + 1)
    ad.X = ad.layers["counts"].copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    # HVGs per sample on raw counts (span 0.5: the default fit failed on this data)
    sc.pp.highly_variable_genes(ad, n_top_genes=args.n_hvgs, flavor="seurat_v3",
                                layer="counts", batch_key="sample", span=0.5)

    # PCA on a scaled copy of the HVGs only (scaling makes the matrix dense)
    hvg = ad[:, ad.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.pp.pca(hvg, n_comps=50)
    ad.obsm["X_pca"] = hvg.obsm["X_pca"]
    ad.uns["pca"] = hvg.uns["pca"]
    del hvg

    # UMAP before integration, kept for before/after plots
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=args.n_pcs, use_rep="X_pca")
    sc.tl.umap(ad, random_state=0)
    ad.obsm["X_umap_uncorrected"] = ad.obsm["X_umap"].copy()

    # Harmony on sample (the array), then neighbors, UMAP and clusters in the corrected space
    ad.obsm["X_pca_harmony"] = run_harmony(ad, "sample")
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=args.n_pcs, use_rep="X_pca_harmony")
    sc.tl.umap(ad, random_state=0)
    for res in args.resolutions:
        sc.tl.leiden(ad, resolution=res, key_added=f"leiden_{res}", flavor="igraph",
                     n_iterations=2, random_state=0)
        print(f"Resolution {res}: {ad.obs[f'leiden_{res}'].nunique()} clusters")

    ad.write_h5ad(args.out)
    print(f"Saved {ad.n_obs} cells to {args.out}")


if __name__ == "__main__":
    main()