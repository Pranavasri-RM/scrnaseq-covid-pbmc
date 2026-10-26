#!/usr/bin/env python3
"""Cell-type composition tests and pseudobulk DESeq2 (COVID vs Healthy, donors as replicates), plus Hallmark enrichment."""

import argparse
import os

import gseapy as gp
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DE_CELL_TYPES = ["CD14 monocyte", "CD16 monocyte", "CD4 T", "CD8 T", "NK", "B cell"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="annotated .h5ad")
    p.add_argument("--outdir", required=True, help="folder for the result tables")
    p.add_argument("--gene-sets", default="MSigDB_Hallmark_2020", help="Enrichr library name or path to a .gmt file")
    return p.parse_args()


def composition(ad):
    """Percent of each donor's cells per cell type, and Mann-Whitney + BH tests per cell type."""
    counts = pd.crosstab(ad.obs["donor"], ad.obs["cell_type"])
    props = counts.div(counts.sum(axis=1), axis=0) * 100
    props["condition"] = ad.obs.groupby("donor", observed=True)["condition"].first()

    is_covid = props["condition"] == "COVID"
    rows = []
    for ct in [c for c in props.columns if c != "condition"]:
        cov, hea = props.loc[is_covid, ct], props.loc[~is_covid, ct]
        p = mannwhitneyu(cov, hea, alternative="two-sided").pvalue
        rows.append({"cell_type": ct, "median_covid": cov.median(), "median_healthy": hea.median(), "p": p})
    comp = pd.DataFrame(rows).set_index("cell_type")
    comp["fdr"] = multipletests(comp["p"], method="fdr_bh")[1]
    return props, comp.sort_values("p").round(3)


def pseudobulk_de(ad, cell_type, min_cells=10, min_donors=3):
    """Pseudobulk DESeq2 (COVID vs Healthy) for one cell type, donors as replicates."""
    sub = ad[ad.obs["cell_type"] == cell_type]
    donors = sub.obs["donor"].astype(str)

    n_cells = donors.value_counts()
    kept = n_cells[n_cells >= min_cells].index

    pb = pd.DataFrame(
        {d: np.asarray(sub[donors == d].layers["counts"].sum(axis=0)).ravel() for d in kept},
        index=sub.var_names,
    ).T
    meta = sub.obs.groupby("donor", observed=True)["condition"].first().astype(str).loc[pb.index].to_frame()

    n_per_group = meta["condition"].value_counts()
    if len(n_per_group) < 2 or n_per_group.min() < min_donors:
        print(f"{cell_type}: skipped, donors per group = {n_per_group.to_dict()}")
        return None

    pb = pb.loc[:, (pb >= 10).sum(axis=0) >= n_per_group.min()]

    dds = DeseqDataSet(counts=pb, metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "COVID", "Healthy"], quiet=True)
    stats.summary()

    out = stats.results_df.copy()
    out["gene"] = out.index.map(ad.var["gene_name"])
    out.attrs["n_donors"] = n_per_group.to_dict()
    return out


def significant(r):
    return r[(r["padj"] < 0.05) & (r["log2FoldChange"].abs() > 1)]


def is_ambient_ig(genes):
    """Antibody genes, treated as ambient RNA from plasmablasts outside B cells."""
    return genes.str.match(r"^IG[HKL][VDJCAGEM]") | (genes == "JCHAIN")


def load_gene_sets(gene_sets):
    """A .gmt path is read into a dict (the same form gseapy uses for a downloaded library); a name is passed through."""
    return gp.read_gmt(gene_sets) if os.path.exists(gene_sets) else gene_sets


def enrichment(de_results, gene_sets):
    """Hallmark over-representation for up and down genes, background = genes tested."""
    gene_sets = load_gene_sets(gene_sets)
    rows = []
    for ct, r in de_results.items():
        r = r[~is_ambient_ig(r["gene"].astype(str))]
        background = r["gene"].dropna().astype(str).tolist()
        s = significant(r)
        for direction, genes in [("up", s[s["log2FoldChange"] > 0]), ("down", s[s["log2FoldChange"] < 0])]:
            if len(genes) < 5:
                continue
            enr = gp.enrichr(gene_list=genes["gene"].astype(str).tolist(), gene_sets=gene_sets,
                             background=background, outdir=None)
            res = enr.results
            if not isinstance(res, pd.DataFrame) or res.empty:   # gseapy returns a list when nothing overlaps
                print(f"{ct} {direction}: no pathway overlaps")
                continue
            top = (res.query("`Adjusted P-value` < 0.05")
                      .sort_values(["Adjusted P-value", "Term"])
                      .head(3))
            for _, t in top.iterrows():
                hit = t["Genes"].split(";")
                rows.append({"cell_type": ct, "direction": direction, "pathway": t["Term"],
                             "n_genes": len(hit), "odds_ratio": round(t["Odds Ratio"], 1),
                             "padj": t["Adjusted P-value"], "example_genes": ", ".join(hit[:6])})
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ad = sc.read_h5ad(args.input)

    props, comp = composition(ad)
    props.to_csv(f"{args.outdir}/proportions_per_donor.csv")
    comp.to_csv(f"{args.outdir}/composition_tests.csv")
    print("Composition tests with FDR < 0.05:", ", ".join(comp.index[comp["fdr"] < 0.05]))

    de_results, summary_rows = {}, []
    for ct in DE_CELL_TYPES:
        r = pseudobulk_de(ad, ct)
        if r is None:
            continue
        de_results[ct] = r
        r.sort_values("padj").to_csv(f"{args.outdir}/pseudobulk_{ct.replace(' ', '_')}.csv")
        s = significant(r)
        up, down = s[s["log2FoldChange"] > 0], s[s["log2FoldChange"] < 0]
        summary_rows.append({
            "cell_type": ct,
            "donors (COVID/Healthy)": f"{r.attrs['n_donors'].get('COVID')}/{r.attrs['n_donors'].get('Healthy')}",
            "genes_tested": len(r),
            "sig_up": len(up),
            "sig_down": len(down),
            "top_up": ", ".join(up.nlargest(5, "log2FoldChange")["gene"]),
            "top_down": ", ".join(down.nsmallest(5, "log2FoldChange")["gene"]),
        })
        print(f"{ct}: {len(up)} up, {len(down)} down")

    pd.DataFrame(summary_rows).set_index("cell_type").to_csv(f"{args.outdir}/de_summary.csv")
    enrichment(de_results, args.gene_sets).to_csv(f"{args.outdir}/hallmark_enrichment.csv", index=False)


if __name__ == "__main__":
    main() 