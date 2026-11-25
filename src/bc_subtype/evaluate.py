"""Reference runner for verifying Paul's baseline handoff to classification/evaluation."""

import importlib.metadata
import json
import platform
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score
from threadpoolctl import threadpool_limits

from bc_subtype.baselines import gene_panels, make_pipeline
from bc_subtype.data import CLASSES, load_dataset, make_folds, write_json
from bc_subtype.download import sha256


def evaluate(data_dir, output_dir, config_path):
    config = json.loads(Path(config_path).read_text())
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if (output_dir / "run.json").exists() or (output_dir / "metrics.tsv").exists():
        raise FileExistsError("Use a fresh output directory to preserve run provenance")
    X, samples, genes = load_dataset(data_dir)
    y = samples["subtype"].to_numpy()
    folds = make_folds(samples, config["n_splits"], config["n_repeats"], config["seed"])
    folds.to_csv(output_dir / "folds.tsv", sep="\t", index=False)
    panels = gene_panels(genes, tuple(config["random_gene_sizes"]), config["random_seed"])
    write_json(
        output_dir / "panels.json",
        {name: genes.iloc[indices]["gene_id"].tolist() for name, indices in panels.items()},
    )
    write_json(
        output_dir / "run.json",
        {
            "status": "running",
            "config": config,
            "python": platform.python_version(),
            "packages": {
                name: importlib.metadata.version(name)
                for name in [
                    "numpy",
                    "pandas",
                    "scipy",
                    "scikit-learn",
                    "threadpoolctl",
                ]
            },
            "prepared_provenance_sha256": sha256(Path(data_dir) / "provenance.json"),
            "config_sha256": sha256(config_path),
            "implementation_sha256": {
                path.name: sha256(path) for path in sorted(Path(__file__).parent.glob("*.py"))
            },
        },
    )
    metrics, predictions = [], []
    confusion = {}
    for (repeat, fold), group in folds.groupby(["repeat", "fold"], sort=True):
        test = group["row_index"].to_numpy()
        train = np.setdiff1d(np.arange(len(samples)), test)
        assert set(samples.iloc[train]["patient_id"]).isdisjoint(samples.iloc[test]["patient_id"])
        for panel, columns in panels.items():
            for classifier in config["classifiers"]:
                start = time.perf_counter()
                model = make_pipeline(classifier, config["seed"], config["n_estimators"])
                with threadpool_limits(limits=config.get("threads", 2)), warnings.catch_warnings():
                    # An unconverged result must not silently become a reported baseline.
                    warnings.simplefilter("error", ConvergenceWarning)
                    model.fit(X[np.ix_(train, columns)], y[train])
                    predicted = model.predict(X[np.ix_(test, columns)])
                per_class = f1_score(
                    y[test], predicted, labels=CLASSES, average=None, zero_division=0
                )
                record = {
                    "repeat": int(repeat),
                    "fold": int(fold),
                    "panel": panel,
                    "classifier": classifier,
                    "input_genes": len(columns),
                    "retained_genes": int(model.named_steps["variance"].get_support().sum()),
                    "train_samples": len(train),
                    "test_samples": len(test),
                    "macro_f1": float(
                        f1_score(
                            y[test], predicted, labels=CLASSES, average="macro", zero_division=0
                        )
                    ),
                    "balanced_accuracy": float(balanced_accuracy_score(y[test], predicted)),
                    "seconds": time.perf_counter() - start,
                    **{f"f1_{label}": float(score) for label, score in zip(CLASSES, per_class)},
                }
                metrics.append(record)
                for index, prediction in zip(test, predicted):
                    predictions.append(
                        {
                            "repeat": int(repeat),
                            "fold": int(fold),
                            "panel": panel,
                            "classifier": classifier,
                            "sample_id": samples.iloc[index]["sample_id"],
                            "true_subtype": y[index],
                            "predicted_subtype": prediction,
                        }
                    )
                key = f"{panel}/{classifier}"
                counts = confusion_matrix(y[test], predicted, labels=CLASSES)
                confusion[key] = confusion.get(key, np.zeros((4, 4), dtype=int)) + counts
                pd.DataFrame(metrics).to_csv(output_dir / "metrics.tsv", sep="\t", index=False)
                print(
                    f"repeat={repeat} fold={fold} {panel}/{classifier} "
                    f"macro-F1={record['macro_f1']:.4f} ({record['seconds']:.1f}s)",
                    flush=True,
                )
    scores = pd.DataFrame(metrics)
    summary = (
        scores.groupby(["panel", "classifier"], sort=False)
        .agg(
            splits=("macro_f1", "count"),
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            balanced_accuracy_mean=("balanced_accuracy", "mean"),
            balanced_accuracy_std=("balanced_accuracy", "std"),
        )
        .reset_index()
    )
    reference = summary.loc[summary["panel"] == "full"].set_index("classifier")["macro_f1_mean"]
    summary["macro_f1_delta_from_full"] = summary["macro_f1_mean"] - summary["classifier"].map(
        reference
    )
    summary["within_0_02_of_full"] = summary["macro_f1_delta_from_full"] >= -0.02
    summary.to_csv(output_dir / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(predictions).to_csv(output_dir / "predictions.tsv", sep="\t", index=False)
    write_json(
        output_dir / "confusion_matrices.json",
        {
            "classes": CLASSES,
            "note": "Sum over splits; repeated CV counts each sample per repeat",
            "matrices": {key: value.tolist() for key, value in confusion.items()},
        },
    )
    run = json.loads((output_dir / "run.json").read_text())
    run["status"] = "complete"
    run["fits"] = len(metrics)
    write_json(output_dir / "run.json", run)
    return summary
