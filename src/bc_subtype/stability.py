"""Stability and PAM50-overlap analysis of fold-wise gene rankings.

Works only from a rankings table (as written by ``bc-subtype rank``), so it needs
no expression data. Each (repeat, fold) training split contributes one selected
set per method and subset size k: the genes with rank <= k.
"""

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from bc_subtype.data import panel_reference, write_json

SPLIT = ["repeat", "fold"]


def load_rankings(path):
    return pd.read_csv(path, sep="\t", dtype={"gene_id": str})


def selected_sets(rankings, method, k):
    """Return {(repeat, fold): set of gene IDs ranked <= k} for one method."""
    chosen = rankings.loc[(rankings["method"] == method) & (rankings["rank"] <= k)]
    sets = {key: set(group["gene_id"]) for key, group in chosen.groupby(SPLIT, sort=True)}
    short = [key for key, genes in sets.items() if len(genes) != k]
    if not sets or short:
        raise ValueError(f"{method}: splits without exactly {k} ranked genes: {short or 'all'}")
    return sets


def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else 1.0


def mean_pairwise_jaccard(sets):
    """Mean Jaccard similarity over all pairs of selected sets (1 = identical)."""
    sets = list(sets)
    if len(sets) < 2:
        raise ValueError("Need at least two selected sets")
    return float(np.mean([jaccard(a, b) for a, b in itertools.combinations(sets, 2)]))


def overlap_pvalue(overlap, k, n_universe, n_panel):
    """P(X >= overlap) when k genes are drawn at random from n_universe genes."""
    return float(hypergeom.sf(overlap - 1, n_universe, n_panel, k))


def selection_frequency(rankings, method, k):
    """How many training splits selected each gene among its top k."""
    sets = selected_sets(rankings, method, k)
    counts = pd.Series([g for genes in sets.values() for g in genes]).value_counts()
    symbols = rankings.drop_duplicates("gene_id").set_index("gene_id")["symbol"]
    pam50 = set(panel_reference()["gene_id"])
    table = pd.DataFrame(
        {
            "method": method,
            "k": k,
            "gene_id": counts.index,
            "symbol": symbols.reindex(counts.index).to_numpy(),
            "in_pam50": counts.index.isin(pam50),
            "splits_selected": counts.to_numpy(),
            "frequency": counts.to_numpy() / len(sets),
        }
    )
    return table.sort_values(["splits_selected", "gene_id"], ascending=[False, True])


def summarize(rankings, n_universe, sizes=(10, 25, 50, 100, 500)):
    """One row per (method, k): stability, PAM50 overlap and agreement across methods."""
    pam50 = set(panel_reference()["gene_id"])
    methods = sorted(rankings["method"].unique())
    rows = []
    for k in sizes:
        by_method = {m: selected_sets(rankings, m, k) for m in methods}
        for method, sets in by_method.items():
            overlaps = np.array([len(genes & pam50) for genes in sets.values()])
            row = {
                "method": method,
                "k": k,
                "splits": len(sets),
                "mean_pairwise_jaccard": mean_pairwise_jaccard(sets.values()),
                "genes_in_all_splits": len(set.intersection(*sets.values())),
                "genes_in_any_split": len(set.union(*sets.values())),
                "pam50_overlap_mean": overlaps.mean(),
                "pam50_overlap_min": int(overlaps.min()),
                "pam50_overlap_max": int(overlaps.max()),
                "pam50_expected_random": k * len(pam50) / n_universe,
                "pam50_overlap_median_p": float(
                    np.median([overlap_pvalue(o, k, n_universe, len(pam50)) for o in overlaps])
                ),
            }
            for other, other_sets in by_method.items():
                if other != method:
                    row[f"mean_overlap_with_{other}"] = float(
                        np.mean([len(sets[s] & other_sets[s]) / k for s in sets])
                    )
            rows.append(row)
    return pd.DataFrame(rows)


def analyze(rankings_path, output_dir, qc_path, sizes=(10, 25, 50, 100, 500)):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rankings = load_rankings(rankings_path)
    n_universe = json.loads(Path(qc_path).read_text())["genes"]
    summary = summarize(rankings, n_universe, sizes)
    frequency = pd.concat(
        [
            selection_frequency(rankings, method, k)
            for method in sorted(rankings["method"].unique())
            for k in sizes
        ],
        ignore_index=True,
    )
    summary.to_csv(output_dir / "stability-summary.tsv", sep="\t", index=False)
    frequency.to_csv(output_dir / "selection-frequency.tsv", sep="\t", index=False)
    write_json(
        output_dir / "stability-run.json",
        {"rankings": str(rankings_path), "n_universe": n_universe, "sizes": list(sizes)},
    )
    return summary
