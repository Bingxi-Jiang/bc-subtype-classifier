import numpy as np
import pandas as pd
import pytest

from bc_subtype.data import CLASSES, clean_expression, make_folds, match_samples


def test_matching_excludes_wrong_tissue_normal_missing_and_duplicate_patient():
    samples = pd.DataFrame(
        [
            ["TCGA-AA-0001", "TCGA-AA-0001-01", "Primary"],
            ["TCGA-AA-0001", "TCGA-AA-0001-02", "Recurrent"],
            ["TCGA-AA-0002", "TCGA-AA-0002-01", "Primary"],
            ["TCGA-AA-0003", "TCGA-AA-0003-01", "Primary"],
        ],
        columns=["PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE"],
    )
    patients = pd.DataFrame(
        [
            ["TCGA-AA-0001", "BRCA_LumA"],
            ["TCGA-AA-0002", "BRCA_Normal"],
            ["TCGA-AA-0003", None],
        ],
        columns=["PATIENT_ID", "SUBTYPE"],
    )
    expression_ids = [
        "TCGA-AA-0001-01B-00R-0000-00",
        "TCGA-AA-0001-01A-00R-0000-00",
        "TCGA-AA-0001-02",
        "TCGA-AA-0002-01",
        "TCGA-AA-0003-01",
        "TCGA-AA-9999-01",
        "invalid",
    ]
    metadata, excluded = match_samples(expression_ids, samples, patients)
    assert metadata["sample_id"].tolist() == ["TCGA-AA-0001-01A-00R-0000-00"]
    assert metadata["subtype"].tolist() == ["Luminal A"]
    assert set(excluded["reason"]) == {
        "duplicate_patient",
        "not_primary_tumor",
        "normal_like_label",
        "missing_subtype",
        "unmatched_or_ambiguous_clinical_sample",
        "invalid_barcode",
    }


def test_unknown_subtype_and_mismatched_patient_fail_closed():
    samples = pd.DataFrame(
        {
            "PATIENT_ID": ["TCGA-AA-0001"],
            "SAMPLE_ID": ["TCGA-AA-0001-01"],
            "SAMPLE_TYPE": ["Primary"],
        }
    )
    patients = pd.DataFrame({"PATIENT_ID": ["TCGA-AA-0001"], "SUBTYPE": ["NEW_LABEL"]})
    with pytest.raises(ValueError, match="Unknown subtype"):
        match_samples(["TCGA-AA-0001-01"], samples, patients)
    samples["PATIENT_ID"] = "TCGA-AA-0002"
    with pytest.raises(ValueError, match="mismatch"):
        match_samples(["TCGA-AA-0001-01"], samples, patients)


def test_expression_collapses_by_entrez_before_log_and_keeps_missing():
    frame = pd.DataFrame(
        {
            "Hugo_Symbol": ["OLD", "NEW", "B", "unknown", "invalid"],
            "Entrez_Gene_Id": [100, 100, 200, np.nan, -1],
            "s1": [3, 5, 0, 10, 10],
            "s2": [7, 9, np.nan, 10, 10],
        }
    )
    X, genes, qc = clean_expression(frame, ["s2", "s1"])
    assert genes["gene_id"].tolist() == ["100", "200"]
    np.testing.assert_allclose(X[:, 0], np.log2([9, 5]), rtol=1e-6)
    assert np.isnan(X[0, 1])
    assert X[1, 1] == 0
    assert qc["invalid_entrez_rows_removed"] == 2
    assert qc["duplicate_entrez_rows_collapsed"] == 1
    assert genes.iloc[0]["symbol"] == "NEW|OLD"


def test_negative_expression_is_not_silently_clipped():
    frame = pd.DataFrame({"Hugo_Symbol": ["A"], "Entrez_Gene_Id": [100], "s1": [-0.1]})
    with pytest.raises(ValueError, match="Negative"):
        clean_expression(frame, ["s1"])


def test_cv_is_deterministic_disjoint_stratified_and_covers_each_patient_per_repeat():
    metadata = pd.DataFrame(
        [
            {"sample_id": f"s{i}_{j}", "patient_id": f"p{i}_{j}", "subtype": label}
            for i, label in enumerate(CLASSES)
            for j in range(10)
        ]
    )
    folds = make_folds(metadata, n_splits=5, n_repeats=5)
    pd.testing.assert_frame_equal(folds, make_folds(metadata, n_splits=5, n_repeats=5))
    assert len(folds) == 5 * len(metadata)
    for _, repeat in folds.groupby("repeat"):
        assert not repeat["patient_id"].duplicated().any()
        assert set(repeat["row_index"]) == set(range(len(metadata)))
    for _, group in folds.groupby(["repeat", "fold"]):
        assert group["subtype"].value_counts().to_dict() == dict.fromkeys(CLASSES, 2)
        test = set(group["patient_id"])
        train = set(metadata["patient_id"]) - test
        assert train.isdisjoint(test)
    metadata.loc[1, "patient_id"] = metadata.loc[0, "patient_id"]
    with pytest.raises(ValueError, match="Patient overlap"):
        make_folds(metadata)
