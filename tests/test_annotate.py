import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import anndata

from annotate import T_GENES, NK_GENES, check_markers, split_nk_cd8


def make_adata(expr, genes, clusters):
    """Small AnnData with Ensembl-style IDs, gene names and a leiden_0.5 column."""
    ad = anndata.AnnData(X=sp.csr_matrix(np.asarray(expr, dtype=float)))
    ad.var_names = [f"ENSG{i:05d}" for i in range(len(genes))]
    ad.var["gene_name"] = genes
    ad.obs["leiden_0.5"] = pd.Categorical(clusters)
    return ad


def test_check_markers_passes_when_labels_match():
    # Cluster 0 expresses CD14, cluster 1 expresses MS4A1
    ad = make_adata([[5, 0], [5, 0], [0, 5], [0, 5]], ["CD14", "MS4A1"], ["0", "0", "1", "1"])
    labels = pd.DataFrame({"cell_type": ["CD14 monocyte", "B cell"], "marker": ["CD14", "MS4A1"]},
                          index=pd.Index(["0", "1"], name="cluster"))
    check_markers(ad, labels)


def test_check_markers_fails_when_clusters_are_renumbered():
    # The bug from the first HPC run: same cells, cluster numbers swapped
    ad = make_adata([[5, 0], [5, 0], [0, 5], [0, 5]], ["CD14", "MS4A1"], ["1", "1", "0", "0"])
    labels = pd.DataFrame({"cell_type": ["CD14 monocyte", "B cell"], "marker": ["CD14", "MS4A1"]},
                          index=pd.Index(["0", "1"], name="cluster"))
    with pytest.raises(AssertionError, match="peaks in cluster"):
        check_markers(ad, labels)


def test_split_nk_cd8_uses_marker_scores():
    rng = np.random.default_rng(0)
    n_fill = 200                                     # background genes, needed for score_genes' control genes
    genes = T_GENES + NK_GENES + [f"GENE{i}" for i in range(n_fill)]
    n = 40
    expr = np.zeros((n, len(genes)))
    expr[:20, :len(T_GENES)] = 3                     # first 20 cells: T-cell genes on
    expr[20:, len(T_GENES):len(T_GENES) + len(NK_GENES)] = 3   # last 20: NK genes on
    expr[:, -n_fill:] = rng.uniform(0, 2, size=(n, n_fill))

    ad = make_adata(expr, genes, ["2"] * n)
    ad.obs["cell_type"] = "NK / CD8 T"
    ad.obsp["connectivities"] = sp.csr_matrix((n, n))  # no neighbors: each cell keeps its own score

    split_nk_cd8(ad)
    assert (ad.obs["cell_type"].iloc[:20] == "CD8 T").all()
    assert (ad.obs["cell_type"].iloc[20:] == "NK").all()