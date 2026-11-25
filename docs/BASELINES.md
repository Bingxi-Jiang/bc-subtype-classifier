# Baseline definitions and experimental protocol

## Fixed panels

| Name | Definition |
| --- | --- |
| `full` | Every valid, deduplicated Entrez feature in the prepared matrix |
| `pam50` | All 50 reference Entrez genes, in reference panel order |
| `random_10`, `random_25`, `random_50`, `random_100`, `random_500` | First k genes from one seeded permutation of the full ordered gene universe |

The random panel seed defaults to 184. Smaller random panels are nested in larger ones, have exactly k distinct input genes, and may contain PAM50 genes. They are fixed across folds and selected without labels or expression values. Save the ordered gene IDs, not only the seed, when reporting a result; the runner writes `panels.json`. A change in the source gene universe or order changes a seeded panel.

One random seed is a reproducible control, not an estimate of variation over random gene sets. The final analysis should repeat random controls using additional predeclared seeds if it makes claims about random-panel distributions. Fixed panels are chosen before fold-local constant filtering, so `retained_genes` may be smaller than the nominal panel size.

The PAM50-gene baseline trains the same Logistic Regression/Random Forest models on those genes. It is **not** the original PAM50 nearest-centroid predictor, does not use a Prosigna assay, and does not generate the target labels. Published subtype annotations are the target for every baseline.

## Fold-local preprocessing

For every outer training fold, each model fits a new scikit-learn pipeline:

1. Median imputation from training samples only. All-missing training genes are imputed to zero, then removed by the next step.
2. Remove genes with zero training variance.
3. Logistic Regression only: standardize using training means and standard deviations.
4. Fit the classifier on that training fold, then transform and predict its held-out fold.

Log2 transformation and per-sample duplicate-row averaging occur before splitting because they are fixed, sample-wise operations. No learned selector belongs in global preparation. Kary's ANOVA or L1 selector must be inserted in a pipeline and fitted on training rows. If tuning hyperparameters or selection size, tune in an **inner** CV inside each outer training fold.

## Reference models

| Model | Fixed settings |
| --- | --- |
| Logistic Regression | `C=1`, `lbfgs`, `max_iter=2000`, balanced class weights, standardized features |
| Random Forest | 200 trees (50 for smoke), sqrt feature sampling, balanced class weights, 2 workers |

Both use seed 184. Thread pools are limited to 2 by default. Convergence warnings are errors: a partial run cannot quietly report an unconverged model as successful. These fixed settings establish usable controls. They are not tuned/optimized model claims. Alan can replace or extend the runner while preserving common splits and fold-local fitting.

## Evaluation and artifacts

`configs/baselines.json` declares repeated stratified 5-fold CV with 5 repeats (25 held-out splits). Every baseline and classifier uses identical assignments for paired comparisons. Since exactly one sample per patient is retained, patient leakage is prevented without a separate grouped splitter.

Per split, `metrics.tsv` records macro-F1, balanced accuracy, each class's F1, sample counts, nominal/retained gene counts, and runtime. `predictions.tsv` records held-out predictions. `summary.tsv` reports mean and sample standard deviation over splits, and the difference from the full-gene mean for the **same classifier**. `confusion_matrices.json` sums counts in fixed class order. In repeated CV each patient contributes once per repeat; totals are not unique patient counts.

The proposal's comparability threshold is a macro-F1 decrease of at most 0.02 from the matching full-gene classifier. `within_0_02_of_full` means `panel_mean - full_mean >= -0.02`; a better panel also qualifies. This is a descriptive project rule, not a noninferiority test or evidence of statistical equivalence. Fold scores share training data; their standard deviation is not an independent-sample confidence interval.

`run.json` records completion state, resolved config, package/Python versions, and checksums of config and prepared provenance. Metrics are checkpointed after each fit. A failure leaves a run marked `running`; only a successful terminal pass changes it to `complete`. Use a new output directory for a rerun. No model pickle is needed to reproduce validation; nothing loads untrusted serialized estimators.

The 2-fold smoke configuration is a real-data integration check. Final project comparisons require the proposal config, Kary's selectors/stability analysis, Alan's reviewed evaluation, and agreed reporting. Do not optimize choices on smoke results and then describe them as prespecified.

## Interpretation limits

PAM50 annotations originate from gene expression, so agreement with them measures reproduction of that molecular classification. It is not independent confirmation of disease, treatment response, survival, or clinical subtype accuracy. Both source batch normalization and subtype annotation construction precede our split. External cohort validation, clinical utility, and biological causal claims are outside this setup. Discuss gene overlap and known pathway roles as interpretation hypotheses, with sources and stability evidence, rather than clinical recommendations.
