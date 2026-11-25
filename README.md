# Interpretable Breast Cancer Subtype Classification

CS 184A project: reproduce TCGA-BRCA molecular subtype annotations using minimal gene-expression features. This initial implementation covers **Paul Jiang's Data Processing & Baselines** responsibilities: pinned dataset acquisition, auditable cohort construction, reproducible baseline panels, and a reference integration runner.

The four targets are Luminal A, Luminal B, Basal-like, and HER2-enriched. The input is batch-normalized RNASeqV2 RSEM expression from the [TCGA PanCancer Atlas breast study](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018). See [the data contract](docs/DATA.md) for the exact release, filtering rules, and limitations.

## Quick start

Python 3.11 or newer is required. Keep the clone and environment on a drive with at least 2 GB free. The three source files total about 155 MB; installation also downloads scientific Python packages.

### Windows PowerShell

```powershell
git clone https://github.com/Bingxi-Jiang/bc-subtype-classifier.git
cd bc-subtype-classifier
# Use an installed Python executable; for the configured target machine:
.\scripts\setup.ps1 -Python D:\Miniconda\python.exe
.\.venv\Scripts\bc-subtype.exe baseline --config configs/smoke.json --output-dir results/smoke
```

The setup script creates `.venv`, installs the package and development tools from PyPI, runs unit tests, downloads the checksum-verified source files, and prepares the cohort. It does not require shell activation. If local policy prevents running a script, use the individual commands below without changing system policy.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --index-url https://pypi.org/simple -e ".[dev]"
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\bc-subtype.exe download
.\.venv\Scripts\bc-subtype.exe prepare
```

### macOS / Linux

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --index-url https://pypi.org/simple -e '.[dev]'
.venv/bin/python -m pytest -q
.venv/bin/bc-subtype download
.venv/bin/bc-subtype prepare
.venv/bin/bc-subtype baseline --config configs/smoke.json --output-dir results/smoke
```

Run commands from the repository root. Optional `--raw-dir`, `--out-dir`, and `--data-dir` arguments support other local storage locations. Data, environments, and generated experiment directories are excluded from Git; teammates reconstruct the same cohort with the download/prepare commands.

## What is included

| Component | Output or purpose |
| --- | --- |
| Pinned downloader | Three public DataHub files, verified by byte count and SHA-256 |
| Cohort preparation | One primary tumor per labelled patient; four classes; audited exclusions |
| Gene representation | Entrez IDs, duplicate collapse, `log2(RSEM + 1)` |
| Baselines | Full gene set, all 50 PAM50 genes, nested random panels of 10/25/50/100/500 genes |
| Fold manifest | Identical stratified 5-fold splits with 5 repeats, seed 184 |
| Reference runner | Logistic Regression and Random Forest with fold-local preprocessing |
| Checks | Synthetic edge-case tests plus Linux/Windows CI; no network data download in CI |

The default proposal-sized baseline run uses 350 fits (7 panels x 2 classifiers x 25 splits):

```powershell
.\.venv\Scripts\bc-subtype.exe baseline --config configs/baselines.json --output-dir results/baselines
```

Use a fresh output directory for each run. The smoke configuration uses real data and all panels, but only 2 folds and 50-tree forests. Its purpose is verifying integration, **not reporting final project performance**. See [validation evidence](docs/VALIDATION.md) for what has actually been run.

An exact dependency snapshot for the validated Windows/Python 3.14 environment is in `requirements-lock-windows-py314.txt`. To recreate that environment, install the lock file into a fresh Python 3.14 virtual environment, then run `python -m pip install --no-deps -e .`. On other Python versions/platforms, use the compatible dependency ranges in `pyproject.toml`; the Windows snapshot is not a universal lock.

## Team handoff

- [Data source, schema, and preprocessing](docs/DATA.md)
- [Baseline definitions and evaluation protocol](docs/BASELINES.md)
- [API examples and division of responsibilities](docs/HANDOFF.md)
- [Reproduction and validation record](docs/VALIDATION.md)
- [Scientific sources and data terms](docs/REFERENCES.md)

Kary owns ANOVA/L1 feature selection and gene stability analysis. Alan owns the final classification/evaluation workflow and visualizations. Their modules can consume `load_dataset`, reuse the committed fold design, and insert fitted selectors inside the provided pipelines. The current reference runner is an integration aid with fixed hyperparameters; it does not claim to complete their research work.

This research predicts existing expression-derived subtype annotations. It does not establish independent diagnostic accuracy, clinical utility, or a replacement for the published PAM50 centroid classifier.
