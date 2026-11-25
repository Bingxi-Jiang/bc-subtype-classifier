"""Match patient labels to primary tumors without learning cohort-level statistics."""

import hashlib
import json
import re
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold

from bc_subtype.download import sha256, verify_raw

LABELS = {
    "BRCA_LumA": "Luminal A",
    "BRCA_LumB": "Luminal B",
    "BRCA_Basal": "Basal-like",
    "BRCA_Her2": "HER2-enriched",
}
CLASSES = list(LABELS.values())
BARCODE = re.compile(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-(\d{2})(?:[A-Z].*|$)")


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n"
    )


def panel_reference():
    with files("bc_subtype").joinpath("resources/pam50.tsv").open() as stream:
        return pd.read_csv(stream, sep="\t", dtype=str)


def match_samples(expression_ids, samples, patients):
    """Return one lexicographically first available primary tumor per labelled patient.

    Exact 15-character sample matches take precedence. Full aliquot barcodes may
    resolve by the first 15 characters only when the clinical mapping is unique.
    Patient-only matching is never used to attach an expression profile.
    """
    if patients["PATIENT_ID"].duplicated().any():
        raise ValueError("Duplicate patient records; resolve labels before matching")
    if samples["SAMPLE_ID"].duplicated().any():
        raise ValueError("Duplicate clinical sample identifiers")
    if len(expression_ids) != len(set(expression_ids)):
        raise ValueError("Duplicate expression sample identifiers")
    labels = patients.set_index("PATIENT_ID")["SUBTYPE"].to_dict()
    exact = samples.set_index("SAMPLE_ID")
    prefixes = {}
    for sample in exact.index:
        prefixes.setdefault(sample[:15], []).append(sample)
    records, exclusions = [], []
    for expression_id in sorted(expression_ids):
        match = BARCODE.fullmatch(expression_id)
        reason = None
        clinical_id = expression_id
        if not match:
            reason = "invalid_barcode"
        elif match.group(2) != "01":
            reason = "not_primary_tumor"
        elif clinical_id not in exact.index:
            candidates = prefixes.get(expression_id[:15], [])
            if len(candidates) == 1:
                clinical_id = candidates[0]
            else:
                reason = "unmatched_or_ambiguous_clinical_sample"
        if reason is None:
            clinical = exact.loc[clinical_id]
            patient_id = match.group(1)
            subtype = labels.get(patient_id)
            if clinical["PATIENT_ID"] != patient_id:
                raise ValueError(f"Clinical patient/barcode mismatch: {expression_id}")
            if str(clinical["SAMPLE_TYPE"]).strip().lower() not in {"primary", "primary tumor"}:
                reason = "clinical_sample_not_primary"
            elif subtype == "BRCA_Normal":
                reason = "normal_like_label"
            elif pd.isna(subtype) or subtype in {None, "", "NA"}:
                reason = "missing_subtype"
            elif subtype not in LABELS:
                raise ValueError(f"Unknown subtype {subtype!r}; update explicit mapping")
            else:
                records.append(
                    {
                        "sample_id": expression_id,
                        "clinical_sample_id": clinical_id,
                        "patient_id": patient_id,
                        "source_subtype": subtype,
                        "subtype": LABELS[subtype],
                    }
                )
        if reason:
            exclusions.append({"sample_id": expression_id, "reason": reason})
    metadata = pd.DataFrame(records)
    if metadata.empty:
        raise ValueError("No eligible expression samples")
    duplicate = metadata["patient_id"].duplicated(keep="first")
    exclusions.extend(
        {"sample_id": s, "reason": "duplicate_patient"}
        for s in metadata.loc[duplicate, "sample_id"]
    )
    metadata = metadata.loc[~duplicate].reset_index(drop=True)
    return metadata, pd.DataFrame(exclusions, columns=["sample_id", "reason"])


def clean_expression(expression, sample_ids):
    """Keep unambiguous Entrez genes; fixed duplicate collapse precedes log2(x+1).

    This does not fit an imputer, filter by variance, or scale using the cohort.
    Missing/infinite expression is left missing for training-fold imputation.
    """
    ids = pd.to_numeric(expression["Entrez_Gene_Id"], errors="coerce")
    valid = ids.notna() & (ids > 0) & (ids % 1 == 0)
    frame = expression.loc[valid].copy()
    frame["gene_id"] = ids.loc[valid].astype("int64").astype(str)
    values = frame[sample_ids].apply(pd.to_numeric, errors="raise")
    values = values.replace([np.inf, -np.inf], np.nan)
    if (values < 0).any().any():
        raise ValueError("Negative expression in nonnegative RSEM profile; check source/units")
    values.index = frame["gene_id"].to_numpy()
    # Multiple expression rows for the same Entrez ID are averaged per sample,
    # using a fixed rule independent of all subtype labels and other patients.
    collapsed = values.groupby(level=0, sort=True).mean()
    symbols = frame.groupby("gene_id")["Hugo_Symbol"].agg(
        lambda x: "|".join(sorted({str(s) for s in x if pd.notna(s)}))
    )
    genes = pd.DataFrame(
        {"gene_id": collapsed.index, "symbol": symbols.reindex(collapsed.index).to_numpy()}
    )
    matrix = np.log2(collapsed.to_numpy(dtype=np.float32).T + 1)
    return (
        matrix,
        genes,
        {
            "source_gene_rows": int(len(expression)),
            "invalid_entrez_rows_removed": int((~valid).sum()),
            "duplicate_entrez_rows_collapsed": int(len(frame) - len(collapsed)),
            "genes": int(len(genes)),
            "missing_expression_values": int(np.isnan(matrix).sum()),
            "transformation": "log2(batch-normalized RSEM + 1), float32",
        },
    )


def make_folds(metadata, n_splits=5, n_repeats=5, seed=184):
    if metadata["patient_id"].duplicated().any():
        raise ValueError("Patient overlap: exactly one sample per patient is required")
    if n_splits < 2 or n_repeats < 1:
        raise ValueError("Require n_splits >= 2 and n_repeats >= 1")
    if metadata["subtype"].value_counts().min() < n_splits:
        raise ValueError("Each class must have at least n_splits samples")
    splitter = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    assignments = []
    for split, (_, test) in enumerate(splitter.split(metadata, metadata["subtype"])):
        for index in test:
            row = metadata.iloc[index]
            assignments.append(
                {
                    "repeat": split // n_splits,
                    "fold": split % n_splits,
                    "row_index": int(index),
                    "sample_id": row["sample_id"],
                    "patient_id": row["patient_id"],
                    "subtype": row["subtype"],
                }
            )
    return pd.DataFrame(assignments)


def prepare(raw_dir, out_dir):
    manifest = verify_raw(raw_dir)
    raw_dir, out_dir = Path(raw_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(raw_dir / "data_clinical_sample.txt", sep="\t", comment="#")
    patients = pd.read_csv(raw_dir / "data_clinical_patient.txt", sep="\t", comment="#")
    expression = pd.read_csv(raw_dir / "data_mrna_seq_v2_rsem.txt", sep="\t")
    expression_ids = list(expression.columns[2:])
    metadata, exclusions = match_samples(expression_ids, samples, patients)
    matrix, genes, gene_qc = clean_expression(expression, metadata["sample_id"].tolist())
    if set(metadata["subtype"]) != set(CLASSES):
        raise ValueError("Prepared cohort does not contain all four target classes")
    reference = panel_reference()
    panel = reference.merge(genes, on="gene_id", how="left", validate="one_to_one")
    panel["available"] = panel["gene_id"].isin(genes["gene_id"])
    if not panel["available"].all():
        missing = panel.loc[~panel["available"], "panel_symbol"].tolist()
        raise ValueError(f"Incomplete PAM50 panel: {missing}")
    np.savez_compressed(out_dir / "expression.npz", X=matrix)
    metadata.to_csv(out_dir / "samples.tsv", sep="\t", index=False, lineterminator="\n")
    genes.to_csv(out_dir / "genes.tsv", sep="\t", index=False, lineterminator="\n")
    exclusions.to_csv(out_dir / "exclusions.tsv", sep="\t", index=False, lineterminator="\n")
    panel.to_csv(out_dir / "pam50_coverage.tsv", sep="\t", index=False, lineterminator="\n")
    make_folds(metadata).to_csv(out_dir / "folds.tsv", sep="\t", index=False, lineterminator="\n")
    cohort_fingerprint = hashlib.sha256(
        metadata.to_csv(index=False, lineterminator="\n").encode()
        + genes.to_csv(index=False, lineterminator="\n").encode()
    ).hexdigest()
    qc = {
        "study_id": manifest["study_id"],
        "source_revision": manifest["revision"],
        "clinical_patients": int(len(patients)),
        "clinical_samples": int(len(samples)),
        "expression_samples": len(expression_ids),
        "retained_samples": int(len(metadata)),
        "unique_patients": int(metadata["patient_id"].nunique()),
        "class_counts": {k: int(v) for k, v in metadata["subtype"].value_counts().items()},
        "source_subtype_counts": {
            str(k): int(v) for k, v in patients["SUBTYPE"].fillna("MISSING").value_counts().items()
        },
        "exclusions": {k: int(v) for k, v in exclusions["reason"].value_counts().items()},
        "clinical_samples_without_expression": int(
            (~samples["SAMPLE_ID"].isin([s[:15] for s in expression_ids])).sum()
        ),
        "pam50_genes_available": int(panel["available"].sum()),
        "cohort_fingerprint": cohort_fingerprint,
        **gene_qc,
    }
    write_json(out_dir / "qc.json", qc)
    write_json(
        out_dir / "provenance.json",
        {
            "source": manifest,
            "outputs": {
                name: sha256(out_dir / name)
                for name in [
                    "expression.npz",
                    "samples.tsv",
                    "genes.tsv",
                    "exclusions.tsv",
                    "pam50_coverage.tsv",
                    "folds.tsv",
                    "qc.json",
                ]
            },
        },
    )
    return qc


def load_dataset(data_dir):
    data_dir = Path(data_dir)
    provenance = json.loads((data_dir / "provenance.json").read_text())
    for name, expected in provenance["outputs"].items():
        if sha256(data_dir / name) != expected:
            raise ValueError(f"Prepared artifact changed: {name}; prepare again")
    with np.load(data_dir / "expression.npz", allow_pickle=False) as archive:
        X = archive["X"]
    samples = pd.read_csv(data_dir / "samples.tsv", sep="\t")
    genes = pd.read_csv(data_dir / "genes.tsv", sep="\t", dtype=str)
    if X.shape != (len(samples), len(genes)):
        raise ValueError("Expression dimensions disagree with sample/gene metadata")
    if samples["sample_id"].duplicated().any() or samples["patient_id"].duplicated().any():
        raise ValueError("Prepared data must contain unique samples and patients")
    return X, samples, genes
