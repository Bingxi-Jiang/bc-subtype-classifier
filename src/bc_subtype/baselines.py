"""Fixed full-gene, PAM50-gene and nested random-gene panels plus fold-safe models."""

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from bc_subtype.data import panel_reference


def gene_panels(genes, sizes=(10, 25, 50, 100, 500), seed=184):
    ids = genes["gene_id"].tolist()
    if len(ids) != len(set(ids)):
        raise ValueError("Gene identifiers must be unique")
    positions = {gene: index for index, gene in enumerate(ids)}
    panel_ids = panel_reference()["gene_id"].tolist()
    missing = set(panel_ids) - set(ids)
    if missing:
        raise ValueError(f"Missing PAM50 Entrez genes: {sorted(missing)}")
    panels = {
        "full": np.arange(len(ids)),
        "pam50": np.array([positions[gene] for gene in panel_ids]),
    }
    # Same seed = same permutation; smaller random panels are nested in larger
    # ones. No expression values or subtype labels enter this choice.
    permutation = np.random.default_rng(seed).permutation(len(ids))
    for size in sizes:
        if not isinstance(size, int) or size < 1 or size > len(ids):
            raise ValueError(f"Invalid random panel size {size} for {len(ids)} genes")
        panels[f"random_{size}"] = permutation[:size]
    return panels


def make_pipeline(classifier, seed=184, n_estimators=200):
    steps = [
        ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
        ("variance", VarianceThreshold(threshold=0)),
    ]
    if classifier == "logistic_regression":
        steps += [
            ("scale", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    C=1.0,
                    solver="lbfgs",
                    max_iter=2000,
                    class_weight="balanced",
                    random_state=seed,
                ),
            ),
        ]
    elif classifier == "random_forest":
        steps += [
            (
                "classifier",
                RandomForestClassifier(
                    n_estimators=n_estimators,
                    max_features="sqrt",
                    class_weight="balanced",
                    random_state=seed,
                    n_jobs=2,
                ),
            )
        ]
    else:
        raise ValueError(f"Unknown classifier: {classifier}")
    return Pipeline(steps)
