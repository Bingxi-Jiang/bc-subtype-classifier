# Validation and reproducibility record

Validated on 2026-11-24 (America/Tijuana), Windows x64, Python 3.14.6. This record describes executed work, not planned experiments.

## Actual cohort

The checksum-verified source contains 1,084 clinical patients/samples and 1,082 expression profiles. The 2 clinical samples without expression are not part of the modelling matrix. Among expression profiles, 101 have missing subtype labels and 36 are Normal-like, leaving **945 distinct primary-tumor patients**.

| Target | Patients |
| --- | ---: |
| Luminal A | 499 |
| Luminal B | 197 |
| Basal-like | 171 |
| HER2-enriched | 78 |
| Total | 945 |

The 20,531 source gene rows become **20,505 unique Entrez genes** after collapsing 26 duplicate rows; no invalid-ID rows were removed. There are no missing expression values in the retained matrix. **All 50 PAM50 genes are available**. These measured patient counts agree with the proposal's estimates.

The cohort fingerprint is `6fb6d3077e89edd342aa6fa0fde8aaa3cfef7c85613e9dd444b123e8e3fa2585`. See [QC JSON](validation/qc.json), [prepared-artifact checksums](validation/provenance.json), and [PAM50 coverage](validation/pam50-coverage.tsv).

## Checks completed locally

- `python -m pip check`: no broken requirements.
- `ruff check .` and `ruff format --check .`: pass.
- `pytest -q`: **11 tests pass**, including complete source-to-metrics integration.
- PowerShell setup script parses successfully; its individual installation/data/test commands were executed.
- A distributable wheel built successfully using existing setuptools 82.0.1 / wheel 0.47.0. Both packaged resource files (source manifest and PAM50 reference) are present. A build-isolation attempt encountered a PyPI network timeout; the successful wheel build used the existing local build tools without isolation.
- All three public source files pass SHA-256 and byte-count verification, including the cached-download path.
- Real-data smoke run completed: **28 fits**, 7 panels x 2 classifiers x 2 folds; 50 trees per forest. See [smoke summary](validation/smoke-summary.tsv) and [run metadata](validation/smoke-run.json).
- Proposal-sized reference baseline run completed: **350 fits**, 7 panels x 2 classifiers x 25 splits; 200 trees per forest. See [full summary](validation/baselines-summary.tsv), [per-fold metrics](validation/baselines-metrics.tsv), [run metadata](validation/baselines-run.json), and [aggregated confusion matrices](validation/baselines-confusion-matrices.json).

The committed CI checks Python 3.11/3.14 on Linux/Windows with synthetic data and requires no 155 MB dataset download. Its actual execution state is available on the [GitHub Actions page](https://github.com/Bingxi-Jiang/bc-subtype-classifier/actions); local success alone does not assert that every remote runner passed.

## Completed baseline results

These runs use `configs/baselines.json` with fixed hyperparameters. The entries below are mean macro-F1 and balanced accuracy over **25 paired held-out splits** (5 folds x 5 repeats). LR means Logistic Regression; RF means Random Forest. Standard deviations and class-specific F1 scores are in the linked TSV files.

| Input panel | LR macro-F1 | RF macro-F1 | LR balanced accuracy | RF balanced accuracy |
| --- | ---: | ---: | ---: | ---: |
| Full: 20,505 genes | 0.8932 | 0.8922 | 0.8994 | 0.8896 |
| PAM50: 50 genes | 0.9363 | 0.9204 | 0.9426 | 0.9268 |
| Random: 10 genes | 0.5660 | 0.5774 | 0.6116 | 0.5954 |
| Random: 25 genes | 0.6923 | 0.7247 | 0.7213 | 0.7350 |
| Random: 50 genes | 0.7500 | 0.7812 | 0.7647 | 0.7787 |
| Random: 100 genes | 0.8049 | 0.8166 | 0.8125 | 0.8089 |
| Random: 500 genes | 0.8780 | 0.8651 | 0.8788 | 0.8595 |

Under the proposal's descriptive 0.02 macro-F1 margin, the PAM50-gene baselines qualify for both classifiers, and the 500-gene random panel qualifies for LR. The 500-gene RF control is 0.0271 below its full-gene reference. This comparison is conditional on these fixed settings, one random-panel seed, this cohort, and these paired splits. It does not establish statistical equivalence, the best possible model, or the project's minimum selected gene count; ANOVA/L1 subsets have not yet been evaluated.

A smaller established gene panel outperforming a full-gene fixed model is consistent with reduced irrelevant-feature burden and the expression-derived target definition. The experiment alone does not identify a causal mechanism or prove a diagnostic advantage. The [baseline protocol](BASELINES.md) explains target dependence, upstream normalization, fold dependence, and the absence of external validation.

## Reproduce

```powershell
.\.venv\Scripts\bc-subtype.exe download
.\.venv\Scripts\bc-subtype.exe prepare
.\.venv\Scripts\bc-subtype.exe baseline --config configs/smoke.json --output-dir results/smoke-reproduction
.\.venv\Scripts\bc-subtype.exe baseline --config configs/baselines.json --output-dir results/baselines-reproduction
```

Use the Windows/Python 3.14 dependency lock to match the validated environment. `pyproject.toml` permits compatible versions on other platforms; package changes can alter fold generation, numerical fitting, serialization hashes, or predictions, so compare the recorded runtime metadata and cohort fingerprint. The reference run stores implementation checksums as well as configuration and data-provenance checksums. Runtime duration is hardware dependent; successful fitting took approximately 186 seconds in total on the setup machine.

Only code/configuration, documentation, and aggregate evidence are committed. Full expression matrices, sample/patient IDs, fold manifests, and prediction-level artifacts remain in ignored local directories and are regenerated from the public source. This record does not claim that Kary's selection/stability work, Alan's final visualization/evaluation review, or the team's final report is finished.
