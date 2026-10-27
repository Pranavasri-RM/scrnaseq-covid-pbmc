#!/usr/bin/env python3
"""README figures from the final pipeline results (results/pipeline/)."""

import argparse
import os

import matplotlib
matplotlib.use("Agg") 
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from adjustText import adjust_text
from matplotlib.lines import Line2D

# Same meaning in every figure: blue = healthy, orange = COVID
HEALTHY, COVID = "#2a78d6", "#eb6834"
INK, MUTED, BACKGROUND_DOTS = "#0b0b0b", "#52514e", "#c9c8c2"
HALO = [pe.withStroke(linewidth=3, foreground="white")]   # white outline keeps labels readable over dots

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.dpi": 200, "savefig.bbox": "tight",
})

DE_CELL_TYPES = ["CD14 monocyte", "CD16 monocyte", "CD4 T", "CD8 T", "NK", "B cell"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", default="results/pipeline", help="folder with annotated.h5ad and de/")
    p.add_argument("--outdir", default="figures")
    return p.parse_args()


def condition_legend(ax, **kwargs):
    handles = [Line2D([0], [0], marker="o", ls="", color=HEALTHY, markersize=7, label="Healthy (6 donors)"),
               Line2D([0], [0], marker="o", ls="", color=COVID, markersize=7, label="COVID-19 (6 patients)")]
    ax.legend(handles=handles, frameon=False, **kwargs)


def fig_umap(ad, out):
    """Left: cell types named on the map. Right: the same map colored by condition."""
    xy = ad.obsm["X_umap"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))

    ax = axes[0]
    ax.scatter(xy[:, 0], xy[:, 1], s=1, c=BACKGROUND_DOTS, linewidths=0, rasterized=True)
    texts = []
    for ct, idx in ad.obs.groupby("cell_type", observed=True).indices.items():
        x, y = np.median(xy[idx], axis=0)
        texts.append(ax.text(x, y, ct, fontsize=8, color=INK, ha="center", va="center", path_effects=HALO))
    adjust_text(texts, ax=ax)
    ax.set_title(f"15 cell types ({ad.n_obs:,} cells)")

    ax = axes[1]
    order = np.random.default_rng(0).permutation(ad.n_obs)   # shuffled so neither group hides the other
    colors = np.where(ad.obs["condition"].values == "COVID", COVID, HEALTHY)
    ax.scatter(xy[order, 0], xy[order, 1], s=1, c=colors[order], linewidths=0, rasterized=True)
    ax.set_title("COVID-19 vs healthy")
    condition_legend(ax, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2)

    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    fig.savefig(out)
    plt.close(fig)


def fig_markers(ad, out):
    """Marker dot plot: evidence for each cell-type label."""
    markers = {
        "CD14 monocyte": ["CD14", "LYZ"],
        "CD16 monocyte": ["FCGR3A", "MS4A7"],
        "Dendritic cell": ["FCER1A", "CD1C"],
        "pDC": ["LILRA4"],
        "CD4 T": ["CD3E", "IL7R"],
        "CD8 T": ["CD8A"],
        "NK": ["KLRF1", "NCAM1"],
        "Proliferating": ["MKI67"],
        "B cell": ["MS4A1"],
        "Plasmablast": ["JCHAIN", "MZB1"],
        "Platelet": ["PPBP"],
        "Red blood cell": ["HBB"],
        "Progenitor": ["CD34"],
        "Neutrophil": ["FCGR3B", "CXCR2"],
        "Developing neutrophil": ["CEACAM8"],
    }
    # var_names are Ensembl IDs; gene symbols live in a column of ad.var
    candidates = ["gene_symbols", "gene_symbol", "gene_name", "gene_names", "symbol"]
    sym_col = next((c for c in candidates if c in ad.var.columns), None)
    assert sym_col is not None, f"no gene-symbol column in ad.var (columns: {list(ad.var.columns)})"

    genes = [g for gs in markers.values() for g in gs]   # flat list: no group brackets above the plot
    dp = sc.pl.dotplot(ad, genes, groupby="cell_type", categories_order=list(markers),
                       gene_symbols=sym_col, cmap="Blues", standard_scale="var", return_fig=True,
                       colorbar_title="Mean expression\n(scaled 0–1 per gene)",
                       size_title="Fraction of cells\nexpressing (%)")
    dp.make_figure()                                   # build the plot so its axes can be edited
    main_ax = dp.get_axes()["mainplot_ax"]             # the axes that hold the dots
    main_ax.set_title("Marker genes per cell type (pipeline labels)", pad=12)
    dp.fig.savefig(out)
    plt.close("all")


def fig_composition(props, comp, out):
    """One dot per donor for the cell types that shift most; line = median."""
    types = ["Plasmablast", "Platelet", "Developing neutrophil", "CD14 monocyte", "pDC", "NK"]
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5))
    for ax, ct in zip(axes.flat, types):
        for i, (cond, color) in enumerate([("Healthy", HEALTHY), ("COVID", COVID)]):
            v = props.loc[props["condition"] == cond, ct].values
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=40, color=color,
                       edgecolor="white", linewidth=0.8, zorder=3)
            ax.hlines(np.median(v), i - 0.25, i + 0.25, color=INK, lw=2, zorder=4)
        ax.set_xticks([0, 1], ["Healthy", "COVID-19"])
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylim(bottom=0)
        ax.set_title(f"{ct}\nFDR = {comp.loc[ct, 'fdr']:.3f}", fontsize=10)
        ax.set_ylabel("% of donor's cells")
    fig.suptitle("Cell-type proportions per donor (Mann-Whitney U, Benjamini-Hochberg across 15 cell types)", y=1.0)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_volcano(res, out):
    """CD14 monocytes, pseudobulk DESeq2; orange = higher in COVID, blue = higher in healthy."""
    res = res.dropna(subset=["padj"])
    lfc = res["log2FoldChange"]
    y = -np.log10(res["padj"].clip(lower=1e-300))
    sig = (res["padj"] < 0.05) & (lfc.abs() > 1)
    up, down = sig & (lfc > 0), sig & (lfc < 0)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(lfc[~sig], y[~sig], s=4, c=BACKGROUND_DOTS, linewidths=0, rasterized=True)
    ax.scatter(lfc[up], y[up], s=9, c=COVID, linewidths=0, label=f"Higher in COVID-19 ({up.sum()})")
    ax.scatter(lfc[down], y[down], s=9, c=HEALTHY, linewidths=0, label=f"Higher in healthy ({down.sum()})")
    ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c=MUTED)
    for v in (-1, 1):
        ax.axvline(v, ls="--", lw=0.8, c=MUTED)

    key_genes = ["IFI44L", "IFITM3", "MX1", "ISG15",            # interferon response
                 "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74",     # antigen presentation (MHC class II)
                 "IL1B", "TNF"]                                 # inflammatory cytokines
    texts = []
    for g in key_genes:
        hit = res[res["gene"] == g]
        if len(hit):
            r = hit.iloc[0]
            gx, gy = r["log2FoldChange"], -np.log10(r["padj"])
            ax.scatter(gx, gy, s=30, facecolors="none", edgecolors=INK, linewidths=1, zorder=4)
            texts.append(ax.text(gx, gy, g, fontsize=8, color=INK, path_effects=HALO))
    adjust_text(texts, ax=ax, expand=(1.4, 1.8), arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))

    ax.set_xlabel("log2 fold change (COVID-19 / healthy)")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title("CD14 monocytes: pseudobulk DESeq2, 6 vs 6 donors")
    ax.legend(frameon=False, loc="upper left")
    fig.savefig(out)
    plt.close(fig)


def fig_enrichment(enr, out):
    """Hallmark pathways per cell type: dot size = genes in the overlap, color = significance."""
    enr = enr[enr["n_genes"] >= 5].copy()   # pathways driven by only a few genes left out of the figure
    enr["label"] = enr["pathway"] + np.where(enr["direction"] == "down", " (down)", "")
    rows = enr.groupby("label")["padj"].min().sort_values().index.tolist()   # most significant first
    x = enr["cell_type"].map({c: i for i, c in enumerate(DE_CELL_TYPES)})
    y = enr["label"].map({r: i for i, r in enumerate(rows)})

    fig, ax = plt.subplots(figsize=(8, 0.45 * len(rows) + 2))
    pts = ax.scatter(x, y, s=enr["n_genes"] * 6, c=-np.log10(enr["padj"]), cmap="Blues", vmin=0,
                     edgecolor=MUTED, linewidth=0.5)
    ax.set_xticks(range(len(DE_CELL_TYPES)), DE_CELL_TYPES, rotation=30, ha="right")
    ax.set_yticks(range(len(rows)), rows)
    ax.set_xlim(-0.6, len(DE_CELL_TYPES) - 0.4)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.grid(color="#e6e5e0", lw=0.6)
    ax.set_axisbelow(True)
    fig.colorbar(pts, ax=ax, label="-log10 adjusted p", shrink=0.7)
    for n in (10, 25, 50):
        ax.scatter([], [], s=n * 6, c="white", edgecolor=MUTED, label=f"{n} genes")
    ax.legend(title="Genes in overlap", frameon=False, loc="upper left", bbox_to_anchor=(1.3, 1), labelspacing=1.5)
    ax.set_title("Hallmark pathways enriched among genes changed in COVID-19\n"
                 "(up-regulated unless marked 'down'; top 3 per cell type, at least 5 genes)")
    fig.savefig(out)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    de = f"{args.results}/de"

    ad = sc.read_h5ad(f"{args.results}/annotated.h5ad")
    fig_umap(ad, f"{args.outdir}/umap_cell_types.png")
    fig_markers(ad, f"{args.outdir}/marker_dotplot.png")

    props = pd.read_csv(f"{de}/proportions_per_donor.csv", index_col=0)
    comp = pd.read_csv(f"{de}/composition_tests.csv", index_col=0)
    fig_composition(props, comp, f"{args.outdir}/composition.png")

    fig_volcano(pd.read_csv(f"{de}/pseudobulk_CD14_monocyte.csv", index_col=0),
                f"{args.outdir}/cd14_volcano.png")
    fig_enrichment(pd.read_csv(f"{de}/hallmark_enrichment.csv"), f"{args.outdir}/hallmark_enrichment.png")
    print(f"Figures written to {args.outdir}/")


if __name__ == "__main__":
    main()