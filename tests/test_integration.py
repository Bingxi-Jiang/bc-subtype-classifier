import json

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import f1_score

from bc_subtype import download
from bc_subtype.data import CLASSES, LABELS, load_dataset, panel_reference, prepare
from bc_subtype.evaluate import evaluate
from bc_subtype.selection import rank_folds


def test_verified_sources_to_aligned_data_and_auditable_held_out_metrics(tmp_path, monkeypatch):
    raw = tmp_path / "raw"
    raw.mkdir()
    patients = pd.DataFrame(
        [
            {"PATIENT_ID": f"TCGA-AA-{i * 10 + j:04}", "SUBTYPE": label}
            for i, label in enumerate(LABELS)
            for j in range(10)
        ]
    )
    samples = patients[["PATIENT_ID"]].copy()
    samples["SAMPLE_ID"] = samples["PATIENT_ID"] + "-01"
    samples["SAMPLE_TYPE"] = "Primary"
    patients.to_csv(raw / "data_clinical_patient.txt", sep="\t", index=False)
    samples.to_csv(raw / "data_clinical_sample.txt", sep="\t", index=False)
    reference = panel_reference()
    expression = reference.rename(
        columns={"panel_symbol": "Hugo_Symbol", "gene_id": "Entrez_Gene_Id"}
    )
    rng = np.random.default_rng(184)
    for index, sample in enumerate(samples["SAMPLE_ID"]):
        values = rng.uniform(0, 1, len(reference))
        values[index // 10] += 20
        expression[sample] = values
    # Reverse the source columns: prepared ordering must come from its metadata.
    expression = expression[list(expression.columns[:2]) + list(expression.columns[2:][::-1])]
    expression.to_csv(raw / "data_mrna_seq_v2_rsem.txt", sep="\t", index=False)
    manifest = {
        "study_id": "test",
        "revision": "synthetic",
        "files": {
            p.name: {"bytes": p.stat().st_size, "sha256": download.sha256(p)} for p in raw.iterdir()
        },
    }
    monkeypatch.setattr(download, "source_manifest", lambda: manifest)
    processed = tmp_path / "processed"
    qc = prepare(raw, processed)
    X, metadata, genes = load_dataset(processed)
    assert X.shape == (40, 50)
    assert qc["class_counts"] == dict.fromkeys(CLASSES, 10)
    first_gene_id = str(reference.iloc[0]["gene_id"])
    column = genes["gene_id"].tolist().index(first_gene_id)
    expected = np.log2(
        expression.set_index("Entrez_Gene_Id").loc[first_gene_id, metadata.iloc[0]["sample_id"]] + 1
    )
    np.testing.assert_allclose(X[0, column], expected, rtol=1e-6)
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "n_splits": 2,
                "n_repeats": 1,
                "seed": 184,
                "random_seed": 184,
                "random_gene_sizes": [10, 25],
                "classifiers": ["logistic_regression"],
                "n_estimators": 5,
                "threads": 1,
            }
        )
    )
    output = tmp_path / "run"
    summary = evaluate(processed, output, config)
    assert len(summary) == 4
    assert (summary["splits"] == 2).all()
    predictions = pd.read_csv(output / "predictions.tsv", sep="\t")
    scores = pd.read_csv(output / "metrics.tsv", sep="\t")
    for (panel, fold), group in predictions.groupby(["panel", "fold"]):
        row = scores.loc[(scores["panel"] == panel) & (scores["fold"] == fold)].iloc[0]
        assert row["macro_f1"] == pytest.approx(
            f1_score(
                group["true_subtype"], group["predicted_subtype"], labels=CLASSES, average="macro"
            )
        )
    for _, group in predictions.groupby("panel"):
        assert len(group) == group["sample_id"].nunique() == 40
    matrices = json.loads((output / "confusion_matrices.json").read_text())["matrices"]
    assert all(np.asarray(matrix).sum() == 40 for matrix in matrices.values())
    assert json.loads((output / "run.json").read_text())["status"] == "complete"
    with pytest.raises(FileExistsError):
        evaluate(processed, output, config)
    # Each class has one gene shifted by +20 (reference rows 0-3); both selectors
    # must recover exactly those genes as the top 4 in every training fold.
    selection = tmp_path / "selection"
    rankings = rank_folds(processed, selection, config, n_rank=4, jobs=1)
    signal = set(reference["gene_id"].astype(str).iloc[:4])
    for _, group in rankings.groupby(["repeat", "fold", "method"]):
        assert set(group["gene_id"].astype(str)) == signal
        assert group["rank"].tolist() == [1, 2, 3, 4]
    assert set(rankings["method"]) == {"anova", "l1"} and rankings["in_pam50"].all()
    assert json.loads((selection / "run.json").read_text())["status"] == "complete"
    with pytest.raises(FileExistsError):
        rank_folds(processed, selection, config)
    (processed / "samples.tsv").write_text("corrupt")
    with pytest.raises(ValueError, match="Prepared artifact changed"):
        load_dataset(processed)
