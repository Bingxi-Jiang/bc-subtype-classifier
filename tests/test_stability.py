import pandas as pd
import pytest
from scipy.stats import hypergeom

from bc_subtype.data import panel_reference
from bc_subtype.stability import (
    mean_pairwise_jaccard,
    overlap_pvalue,
    selected_sets,
    selection_frequency,
    summarize,
)

PAM = panel_reference()["gene_id"].tolist()


def rankings_from(lists):
    """lists: {(method, repeat, fold): [gene ids in rank order]}."""
    rows = [
        {"method": m, "repeat": r, "fold": f, "rank": i + 1, "gene_id": g, "symbol": f"S{g}"}
        for (m, r, f), genes in lists.items()
        for i, g in enumerate(genes)
    ]
    return pd.DataFrame(rows)


def test_jaccard_and_hypergeometric_against_hand_computed_values():
    assert mean_pairwise_jaccard([{1, 2}, {1, 2}]) == 1.0
    assert mean_pairwise_jaccard([{1, 2}, {2, 3}, {4, 5}]) == pytest.approx((1 / 3 + 0 + 0) / 3)
    with pytest.raises(ValueError):
        mean_pairwise_jaccard([{1}])
    assert overlap_pvalue(0, 10, 100, 5) == 1.0
    assert overlap_pvalue(2, 10, 100, 5) == pytest.approx(1 - hypergeom.cdf(1, 100, 5, 10))


def test_sets_frequency_and_summary_on_known_rankings():
    rankings = rankings_from(
        {
            ("anova", 0, 0): [PAM[0], PAM[1], "x1"],
            ("anova", 0, 1): [PAM[0], PAM[1], "x2"],
            ("l1", 0, 0): [PAM[0], "x3", "x1"],
            ("l1", 0, 1): ["x4", PAM[0], "x1"],
        }
    )
    assert selected_sets(rankings, "anova", 2) == {
        (0, 0): {PAM[0], PAM[1]},
        (0, 1): {PAM[0], PAM[1]},
    }
    with pytest.raises(ValueError, match="exactly 4"):
        selected_sets(rankings, "anova", 4)

    freq = selection_frequency(rankings, "l1", 2).set_index("gene_id")
    assert freq.loc[PAM[0], "splits_selected"] == 2 and freq.loc[PAM[0], "in_pam50"]
    assert freq.loc["x3", "frequency"] == 0.5 and not freq.loc["x3", "in_pam50"]

    summary = summarize(rankings, n_universe=1000, sizes=(2,)).set_index("method")
    assert summary.loc["anova", "mean_pairwise_jaccard"] == 1.0
    assert summary.loc["l1", "mean_pairwise_jaccard"] == pytest.approx(1 / 3)
    assert summary.loc["anova", "genes_in_all_splits"] == 2
    assert summary.loc["l1", "genes_in_any_split"] == 3
    assert summary.loc["anova", "pam50_overlap_mean"] == 2.0
    assert summary.loc["l1", "pam50_overlap_mean"] == 1.0
    assert summary.loc["anova", "pam50_expected_random"] == pytest.approx(2 * 50 / 1000)
    # Split (0,0): {P0,P1} vs {P0,x3} -> 1/2; split (0,1): {P0,P1} vs {x4,P0} -> 1/2
    assert summary.loc["anova", "mean_overlap_with_l1"] == 0.5
