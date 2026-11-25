# Team integration handoff

## Responsibilities from the project proposal

| Team member | Ownership | Initial setup status |
| --- | --- | --- |
| Paul Jiang | Acquire/prepare TCGA-BRCA, barcode/label matching, full/PAM50/random baselines; interpretation/report support | Implemented here, with validation evidence |
| Kary Zheng | ANOVA and L1 feature selection; subset-size comparison; selection stability and PAM50 overlap | Extension points and data ready; research implementation remains with Kary |
| Alan He | Logistic Regression/Random Forest, repeated CV, metrics, visualizations | Baseline reference runner supplied for integration; final evaluation remains with Alan |
| All | Experimental design, integration, interpretation, final report | Shared decisions required as the project develops |

Names identify code ownership only. The proposal's student IDs, email addresses, and PDF are not needed to reproduce this repository.

## Consume Paul's data

```python
from bc_subtype.data import load_dataset

X, samples, genes = load_dataset("data/processed")
y = samples["subtype"].to_numpy()
# X[i, j] corresponds to samples.iloc[i] and genes.iloc[j].
# X is already log2(RSEM + 1), but has not been imputed or standardized.
```

```python
from bc_subtype.baselines import gene_panels, make_pipeline

panels = gene_panels(genes, sizes=(10, 25, 50, 100, 500), seed=184)
columns = panels["pam50"]
model = make_pipeline("logistic_regression", seed=184)
# Fit model only with an outer training fold, never with all X before CV.
```

Use gene IDs for saved selections and biological joins. Human-readable symbols can be looked up in `genes.tsv` or `pam50_coverage.tsv`. Store the exact gene universe/cohort provenance with experiments.

## Reuse the paired fold assignments

```python
import numpy as np
import pandas as pd

folds = pd.read_csv("data/processed/folds.tsv", sep="\t")
for (repeat, fold), held_out in folds.groupby(["repeat", "fold"]):
    test = held_out["row_index"].to_numpy()
    train = np.setdiff1d(np.arange(len(samples)), test)
    model = make_pipeline("logistic_regression")
    model.fit(X[train][:, columns], y[train])
    predicted = model.predict(X[test][:, columns])
```

The manifest contains each test assignment once per repeat. Every other patient in that repeat is a training patient. Never randomly split individual predictions or mix fold manifests from different cohort fingerprints.

## Add a training-only selector

For an ANOVA extension, insert `SelectKBest(f_classif, k=...)` before the classifier (after scaling for this pipeline). Keep `full`, `pam50`, and `random_*` runs as controls. The example is an integration sketch, not the project's final ANOVA implementation.

```python
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.pipeline import Pipeline

base = make_pipeline("logistic_regression")
selected = Pipeline(base.steps[:-1] + [("select", SelectKBest(f_classif, k=25))] + base.steps[-1:])
```

Recover original IDs by applying the fitted variance mask to the input panel's genes, then applying the selector's mask. Do not interpret selector offsets directly against the unfiltered full gene table. L1 selection also belongs inside each training fold. Prespecify how to obtain exact subset sizes from coefficient rankings and how to resolve ties; these are Kary's design choices.

For stability analysis, save selected Entrez sets for each repeat/fold and calculate selection frequencies and Jaccard overlap using those sets. Retain per-fold selections rather than only a final model. PAM50 overlap uses the packaged 50-ID list. Interpret consistently selected genes with primary biological literature; no biological enrichment or stability conclusions have been fabricated in the setup.

## Commands for contributors

```powershell
.\.venv\Scripts\python.exe -m ruff check .
.\.venv\Scripts\python.exe -m ruff format --check .
.\.venv\Scripts\python.exe -m pytest -q
```

Commit code/configs/docs and concise aggregate evidence. Keep source data, patient-level artifacts, generated predictions, virtual environments, and large experiment files in the ignored local directories. The downloader provides reproducible access without committing a 155 MB source file. If the source snapshot must change, update its commit/checksums and regenerate QC, fold manifests, panels, and validation evidence together.
