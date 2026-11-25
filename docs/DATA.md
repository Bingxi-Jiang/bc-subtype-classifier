# Data acquisition and preprocessing contract

## Source snapshot

| Property | Value |
| --- | --- |
| Study | `brca_tcga_pan_can_atlas_2018` |
| Provider | cBioPortal DataHub, derived from TCGA PanCancer Atlas |
| DataHub commit | `dda832891555904730395d476c5c888b2e1553e1` |
| Profile | `rna_seq_v2_mrna` |
| Measurement | RSEM, batch normalized from Illumina HiSeq_RNASeqV2 |
| Label field | `SUBTYPE` in patient clinical data |
| Acquisition | Public Git LFS media URLs pinned to the commit; no token required |

The manifest is packaged at `src/bc_subtype/resources/sources.json`. Every raw file is checked against its upstream Git LFS SHA-256 and byte count before processing. Git LFS pointer text, incomplete transfers, and corrupt caches fail validation. Downloads write to a `.part` file and only rename after verification. Re-running download reuses valid files. Move a corrupt file aside before re-downloading; partial files are overwritten on retry.

Only the expression file and the two clinical files are downloaded. The full study archive, mutations, survival outcomes, and unrelated clinical variables are unnecessary for this task. The original proposal PDF is not distributed in this repository.

## Labels and sample selection

| Source label | Model label | Decision |
| --- | --- | --- |
| `BRCA_LumA` | `Luminal A` | Keep |
| `BRCA_LumB` | `Luminal B` | Keep |
| `BRCA_Basal` | `Basal-like` | Keep |
| `BRCA_Her2` | `HER2-enriched` | Keep |
| `BRCA_Normal` | - | Exclude Normal-like |
| Missing | - | Exclude |
| Any new label | - | Stop for explicit review |

1. Start with samples that have expression columns, sorted by expression barcode.
2. Require a valid TCGA sample/aliquot barcode and sample type code `01` (primary solid tumor).
3. Match the expression sample to the exact clinical `SAMPLE_ID`. A full aliquot barcode may resolve to the first 15 characters only if that clinical prefix maps uniquely. Never attach expression using only a patient ID.
4. Require agreement between the barcode-derived patient and the clinical `PATIENT_ID`, and clinical sample type `Primary` or `Primary Tumor`.
5. Join to the unique patient record for `SUBTYPE`, using only the four explicit target mappings.
6. Keep the lexicographically first eligible expression barcode per patient. This fixed choice does not inspect expression, subtype performance, or missingness. Record later samples as duplicate-patient exclusions.

Duplicate clinical sample/patient identifiers, conflicting patient mappings, and unknown labels stop the process. `exclusions.tsv` records an exclusion reason for each expression column not retained. `qc.json` separately reports clinical samples without expression. Expected proposal counts are hypotheses to verify, not constants used to select records. Actual counts appear in [VALIDATION.md](VALIDATION.md) and the committed aggregate QC record.

## Gene processing

The canonical feature key is the positive integer Entrez Gene ID, stored as a string. Gene symbols are annotations only. Rows with absent, noninteger, zero, or negative IDs are excluded. Multiple rows for one Entrez ID are averaged separately for each sample **on the original RSEM scale**, then transformed by `log2(x + 1)`. All available symbols for that ID are retained as a sorted `|`-delimited annotation. This deterministic collapse is label-independent and does not estimate a parameter across patients.

Non-numeric measurements fail parsing. Infinite values become missing. Negative values fail because they contradict this nonnegative RSEM profile. Missing expression remains `NaN` in the saved matrix. No cohort-wide imputation, variance-based filtering, standardization, ANOVA, or feature selection is performed in preparation.

The upstream expression has already undergone normalization/batch processing. That upstream procedure cannot be refitted within our CV folds. This limitation must be acknowledged when interpreting internal validation. These are not raw read counts, TPM values, or z-scores, and this workflow must not be applied unchanged to another measurement profile.

## PAM50 reference

`src/bc_subtype/resources/pam50.tsv` contains 50 unique Entrez IDs and original panel symbols extracted from `pam50$centroids.map` in the [genefu reference dataset](https://github.com/bhklab/genefu/blob/76fc3ca34eee4b55a01219d37da95b7faddaf28f/data/pam50.rda). Its pinned commit is `76fc3ca34eee4b55a01219d37da95b7faddaf28f`. The original R object was read once during setup using Python `rdata`; R and `rdata` are not runtime dependencies.

Joining by Entrez ID resolves old aliases such as `CDCA1`/`NUF2` (83540), `KNTC2`/`NDC80` (10403), and `ORC6L`/`ORC6` (23594). `pam50_coverage.tsv` pairs the original panel symbol with the observed expression annotation. Preparation fails if any of the 50 genes is missing, instead of silently running a smaller panel.

The downloaded reference R object had SHA-256 `766cbeff768d725fd1c800b793d3a0e4d8b71b3fa3d15950410ef076352e7121`.

## Prepared artifacts

All files below live under `data/processed/` by default and are regenerated locally.

| File | Schema / meaning |
| --- | --- |
| `expression.npz` | Key `X`: samples x genes, float32, log2-transformed; may contain NaN |
| `samples.tsv` | Row order of `X`: sample ID, matched clinical sample ID, patient ID, source subtype, model subtype |
| `genes.tsv` | Column order of `X`: `gene_id`, `symbol` |
| `exclusions.tsv` | Expression `sample_id`, exclusion `reason` |
| `pam50_coverage.tsv` | `panel_symbol`, `gene_id`, observed `symbol`, `available` |
| `folds.tsv` | Test-fold assignments: repeat/fold/row index/sample/patient/subtype; indices start at 0 |
| `qc.json` | Counts, missingness, source revision, transformation, cohort fingerprint |
| `provenance.json` | Raw-source manifest and checksums of all prepared artifacts |

For each repeat, every retained patient appears in exactly one test fold. Training rows are the complement of that test fold. The preparation default is 5 folds x 5 repeats with seed 184. The baseline runner emits the assignments for its own configuration so the 2-fold smoke run cannot be confused with the proposal design.

`load_dataset` verifies prepared checksums and dimensions before returning `(X, samples, genes)`. Editing any prepared artifact requires regenerating the cohort and provenance, rather than silently changing its interpretation.
