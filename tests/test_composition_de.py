import numpy as np
import pandas as pd

from composition_de import is_ambient_ig, significant


def test_significant_needs_both_cutoffs():
    r = pd.DataFrame({"padj": [0.01, 0.01, 0.20, np.nan], "log2FoldChange": [2.0, 0.5, 3.0, -4.0]},
                     index=["both", "small_change", "not_significant", "filtered_out"])
    assert list(significant(r).index) == ["both"]


def test_ambient_ig_genes():
    genes = pd.Series(["IGKC", "IGHM", "JCHAIN", "IGLL5", "CD3E", "IGFBP7"])
    assert is_ambient_ig(genes).tolist() == [True, True, True, False, False, False]