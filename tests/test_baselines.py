import numpy as np
import pandas as pd
import pytest

from bc_subtype.baselines import gene_panels, make_pipeline
from bc_subtype.data import panel_reference


def test_pam50_ids_and_random_panels_are_complete_nested_and_reproducible():
    reference = panel_reference()
    assert len(reference) == reference["gene_id"].nunique() == 50
    genes = pd.DataFrame({"gene_id": reference["gene_id"].tolist() + ["999999"]})
    first = gene_panels(genes, sizes=(10, 25, 50))
    second = gene_panels(genes, sizes=(10, 25, 50))
    for key in first:
        np.testing.assert_array_equal(first[key], second[key])
    assert genes.iloc[first["pam50"]]["gene_id"].tolist() == reference["gene_id"].tolist()
    assert set(first["random_10"]) <= set(first["random_25"]) <= set(first["random_50"])
    assert len(first["full"]) == 51
    assert not np.array_equal(
        first["random_10"], gene_panels(genes, sizes=(10,), seed=1)["random_10"]
    )
    with pytest.raises(ValueError, match="Invalid random panel"):
        gene_panels(genes, sizes=(52,))
    with pytest.raises(ValueError, match="Missing PAM50"):
        gene_panels(genes.iloc[1:], sizes=())


@pytest.mark.parametrize("classifier", ["logistic_regression", "random_forest"])
def test_preprocessing_uses_only_training_values_and_handles_empty_constant_genes(classifier):
    X = np.array(
        [
            [0, 10, np.nan, 7],
            [2, 20, np.nan, 7],
            [4, 30, np.nan, 7],
            [6, np.nan, np.nan, 7],
        ]
    )
    model = make_pipeline(classifier, n_estimators=5)
    model.fit(X, ["a", "a", "b", "b"])
    np.testing.assert_allclose(model.named_steps["impute"].statistics_, [3, 20, 0, 7])
    np.testing.assert_array_equal(
        model.named_steps["variance"].get_support(), [True, True, False, False]
    )
    statistics_before = model.named_steps["impute"].statistics_.copy()
    assert len(model.predict([[100000, np.nan, 10000, -500]])) == 1
    np.testing.assert_array_equal(model.named_steps["impute"].statistics_, statistics_before)
    if classifier == "logistic_regression":
        np.testing.assert_allclose(model.named_steps["scale"].mean_, [3, 20])
