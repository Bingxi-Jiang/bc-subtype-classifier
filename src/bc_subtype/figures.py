"""
Reads the per-split metrics written by `bc-subtype baseline ... --rankings ...` and draws:

fig1: macro-F1 versus number of selected genes, one panel per classifier.
fig2: row-normalized confusion matrices for the reference and selected gene panels.
fig3: F1 of each subtype versus number of selected genes.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import NullLocator

SIZES = [10, 25, 50, 100, 500]
MARGIN = 0.02
CLASSIFIERS = [
    ("logistic_regression", "Logistic regression"),
    ("random_forest", "Random forest"),
]
# panel prefix, legend label, color, marker, x shift so error bars at the same size don't overlap
LINES = [
    ("l1", "L1 selection", "#257de0", "o", 0.92),
    ("anova", "ANOVA selection", "#f47230", "s", 1.00),
    ("random", "Random genes", "#8c8c8c", "^", 1.08),
]
PAM50_COLOR = "#00b894"
INK, MUTED, AXIS, GRID, BAND = "#191919", "#666666", "#c2c2c2", "#e7e7e7", "#eeeeee"


def load_metrics(results_dir):
    """One row per CV split, panel and classifier."""
    return pd.read_csv(Path(results_dir) / "metrics.tsv", sep="\t")


def summarize(metrics):
    """Mean and standard deviation of macro-F1 over the CV splits, per panel and classifier."""
    return metrics.groupby(["panel", "classifier"])["macro_f1"].agg(["mean", "std"])


def smallest_within_margin(stats, prefix, classifier):
    """Smallest panel size whose mean macro-F1 is at most MARGIN below the all-gene model."""
    full = stats.loc[("full", classifier), "mean"]
    for size in SIZES:
        if stats.loc[(f"{prefix}_{size}", classifier), "mean"] >= full - MARGIN:
            return size
    return None


def style_axis(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def plot_performance_vs_size(metrics, out_dir):
    stats = summarize(metrics)
    n_genes = int(metrics.loc[metrics["panel"] == "full", "input_genes"].iloc[0])
    first = CLASSIFIERS[0][0]
    n_splits = int(((metrics["panel"] == "full") & (metrics["classifier"] == first)).sum())

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    for ax, (classifier, title) in zip(axes, CLASSIFIERS):
        style_axis(ax)

        # Reference: the all-gene model and the band of panels counted as "comparable"
        full = stats.loc[("full", classifier), "mean"]
        ax.axhspan(full - MARGIN, full, color=BAND, zorder=0, label=f"Within {MARGIN} of all genes")
        ax.axhline(full, color=INK, linewidth=1.2, zorder=1, label=f"All {n_genes:,} genes")

        # One line per selection method or control: mean with ±1 SD over the CV splits
        for prefix, label, color, marker, shift in LINES:
            rows = stats.loc[[(f"{prefix}_{size}", classifier) for size in SIZES]]
            ax.errorbar(
                [size * shift for size in SIZES],
                rows["mean"],
                yerr=rows["std"],
                color=color,
                marker=marker,
                markersize=6,
                linewidth=2,
                elinewidth=1,
                capsize=0,
                label=label,
                zorder=3,
            )

        # PAM50: the 50 genes the labels were defined with
        pam50 = stats.loc[("pam50", classifier)]
        ax.errorbar(
            [50],
            [pam50["mean"]],
            yerr=[pam50["std"]],
            color=PAM50_COLOR,
            marker="*",
            markersize=13,
            linestyle="none",
            elinewidth=1,
            capsize=0,
            label="PAM50 genes (50)",
            zorder=4,
        )

        # Panel title plus the finding computed from the data
        found = []
        for prefix, label, *_ in LINES:
            size = smallest_within_margin(stats, prefix, classifier)
            found.append(f"{label.split()[0]} {size if size else 'not reached'}")
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=INK, pad=20)
        ax.text(
            0,
            1.015,
            f"Smallest panel within {MARGIN} of all genes: " + " · ".join(found),
            transform=ax.transAxes,
            fontsize=9,
            color=MUTED,
            va="bottom",
        )

        ax.set_xscale("log")
        ax.set_xticks(SIZES)
        ax.set_xticklabels([str(size) for size in SIZES])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlim(7.5, 650)
        ax.set_xlabel("Number of genes (log scale)", color=MUTED)

    axes[0].set_ylabel(f"Macro-F1 (mean ± SD over {n_splits} CV splits)", color=MUTED)

    # One shared legend, in reading order
    handles, labels = axes[0].get_legend_handles_labels()
    order = [label for _, label, *_ in LINES] + [
        "PAM50 genes (50)",
        f"All {n_genes:,} genes",
        f"Within {MARGIN} of all genes",
    ]
    by_label = dict(zip(labels, handles))
    fig.legend(
        [by_label[label] for label in order],
        order,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.1),
        ncol=6,
        frameon=False,
        fontsize=9,
    )
    fig.tight_layout()

    path = Path(out_dir) / "fig1_macro_f1_vs_genes.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return path


SHORT_NAMES = {
    "Luminal A": "LumA",
    "Luminal B": "LumB",
    "Basal-like": "Basal",
    "HER2-enriched": "HER2",
}
# panel name in the results, title shown above the matrix
CONFUSION_PANELS = [
    ("full", "All genes"),
    ("pam50", "PAM50 genes (50)"),
    ("l1_25", "L1 top 25"),
    ("anova_25", "ANOVA top 25"),
]


def plot_confusion_matrices(results_dir, metrics, out_dir, classifier="logistic_regression"):
    data = json.loads((Path(results_dir) / "confusion_matrices.json").read_text())
    classes = data["classes"]
    names = [SHORT_NAMES.get(name, name) for name in classes]
    stats = summarize(metrics)

    fig, axes = plt.subplots(1, len(CONFUSION_PANELS), figsize=(14, 4.2), layout="constrained")
    for ax, (panel, title) in zip(axes, CONFUSION_PANELS):
        counts = np.array(data["matrices"][f"{panel}/{classifier}"], dtype=float)
        # Divide each row by its total: row i shows where true class i ended up (row sums to 1)
        rates = counts / counts.sum(axis=1, keepdims=True)
        image = ax.imshow(rates, cmap="Blues", vmin=0, vmax=1)
        for i in range(len(classes)):
            for j in range(len(classes)):
                value = rates[i, j]
                ax.text(
                    j,
                    i,
                    f"{value:.0%}" if value >= 0.1 else f"{value:.1%}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="white" if value > 0.6 else INK,
                )
        ax.set_xticks(range(len(classes)), names)
        ax.set_yticks(range(len(classes)), names)
        ax.tick_params(length=0, colors=MUTED, labelsize=9)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlabel("Predicted subtype", color=MUTED)
        macro_f1 = stats.loc[(panel, classifier), "mean"]
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=INK, pad=18)
        ax.text(
            0,
            1.01,
            f"macro-F1 {macro_f1:.3f}",
            transform=ax.transAxes,
            fontsize=9,
            color=MUTED,
            va="bottom",
        )
    axes[0].set_ylabel("True subtype", color=MUTED)
    colorbar = fig.colorbar(image, ax=axes, fraction=0.02, pad=0.02)
    colorbar.set_label("Share of the true subtype", color=MUTED)
    colorbar.ax.tick_params(colors=MUTED, labelsize=9)
    colorbar.outline.set_visible(False)

    path = Path(out_dir) / f"fig2_confusion_{classifier}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return path


def plot_per_class_f1(results_dir, metrics, out_dir, classifier="logistic_regression"):
    data = json.loads((Path(results_dir) / "confusion_matrices.json").read_text())
    classes = data["classes"]
    # Every patient is tested once per repeat, so a row total divided by the repeats = class size
    counts = np.array(data["matrices"][f"full/{classifier}"])
    class_sizes = counts.sum(axis=1) / metrics["repeat"].nunique()

    rows = metrics[metrics["classifier"] == classifier]
    n_genes = int(rows.loc[rows["panel"] == "full", "input_genes"].iloc[0])
    n_splits = int((rows["panel"] == "full").sum())

    fig, axes = plt.subplots(1, len(classes), figsize=(14, 4.4), sharey=True)
    for ax, name, n_patients in zip(axes, classes, class_sizes):
        style_axis(ax)
        # Same as summarize(), but for this subtype's F1 column instead of macro-F1
        stats = rows.groupby("panel")[f"f1_{name}"].agg(["mean", "std"])

        ax.axhline(
            stats.loc["full", "mean"],
            color=INK,
            linewidth=1.2,
            zorder=1,
            label=f"All {n_genes:,} genes",
        )
        for prefix, label, color, marker, shift in LINES:
            panel_stats = stats.loc[[f"{prefix}_{size}" for size in SIZES]]
            ax.errorbar(
                [size * shift for size in SIZES],
                panel_stats["mean"],
                yerr=panel_stats["std"],
                color=color,
                marker=marker,
                markersize=6,
                linewidth=2,
                elinewidth=1,
                capsize=0,
                label=label,
                zorder=3,
            )
        pam50 = stats.loc["pam50"]
        ax.errorbar(
            [50],
            [pam50["mean"]],
            yerr=[pam50["std"]],
            color=PAM50_COLOR,
            marker="*",
            markersize=13,
            linestyle="none",
            elinewidth=1,
            capsize=0,
            label="PAM50 genes (50)",
            zorder=4,
        )

        ax.set_title(name, loc="left", fontsize=11, fontweight="bold", color=INK, pad=20)
        ax.text(
            0,
            1.015,
            f"{n_patients:.0f} patients",
            transform=ax.transAxes,
            fontsize=9,
            color=MUTED,
            va="bottom",
        )
        ax.set_xscale("log")
        ax.set_xticks(SIZES)
        ax.set_xticklabels([str(size) for size in SIZES])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlim(7.5, 650)
        ax.set_xlabel("Number of genes (log scale)", color=MUTED)

    axes[0].set_ylabel(f"F1 for this subtype (mean ± SD, {n_splits} splits)", color=MUTED)

    handles, labels = axes[0].get_legend_handles_labels()
    order = [label for _, label, *_ in LINES] + ["PAM50 genes (50)", f"All {n_genes:,} genes"]
    by_label = dict(zip(labels, handles))
    fig.legend(
        [by_label[label] for label in order],
        order,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.1),
        ncol=5,
        frameon=False,
        fontsize=9,
    )
    fig.tight_layout()

    path = Path(out_dir) / f"fig3_per_class_f1_{classifier}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    parser = argparse.ArgumentParser(description="Draw evaluation figures from a results directory")
    parser.add_argument("--results", default="results/final", help="Folder with metrics.tsv")
    parser.add_argument("--out", default="figures", help="Folder to save the figures in")
    args = parser.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = load_metrics(args.results)
    print(f"Saved {plot_performance_vs_size(metrics, out_dir)}")
    print(f"Saved {plot_confusion_matrices(args.results, metrics, out_dir)}")
    print(f"Saved {plot_per_class_f1(args.results, metrics, out_dir)}")


if __name__ == "__main__":
    main()
