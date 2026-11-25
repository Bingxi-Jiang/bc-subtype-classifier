"""Training-fold-only gene ranking and top-k selection (Kary's feature selection).

Each ranking function takes the training rows of one fold and returns *every*
column index ordered from most to least informative. A ranking is computed once
per training fold; the top 10/25/50/100/500 genes are then nested prefixes of it.
Ties are always broken by column index so that results are fully reproducible.
"""

import warnings

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import SelectorMixin, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.svm import l1_min_c
from sklearn.utils.validation import check_is_fitted, validate_data


def anova_ranking(X, y):
    """Rank columns by one-way ANOVA F statistic across classes (largest first).

    Returns (order, scores). Columns whose F is undefined (e.g. constant in the
    training fold) receive score -inf and are ranked last.
    """
    # Constant columns are handled explicitly below, so their warnings are noise.
    with np.errstate(divide="ignore", invalid="ignore"), warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Features .* are constant")
        scores, _ = f_classif(X, y)
    scores = np.where(np.isfinite(scores), scores, -np.inf)
    # lexsort uses the LAST key as primary: sort by -score, then by column index.
    order = np.lexsort((np.arange(X.shape[1]), -scores))
    return order, scores


def l1_path_ranking(X, y, n_rank=500, step=1.1, max_steps=300, seed=184):
    """Rank columns by entry order along an L1 logistic-regression path.

    Columns are standardized with this training fold's statistics. Starting at
    the smallest C for which some coefficient can be nonzero, C is multiplied by
    ``step`` until at least ``n_rank`` columns have entered. A column "enters" at
    the first C where any of the one-vs-rest class models gives it a nonzero
    weight. Columns entering at the same C are ordered by their largest absolute
    weight at that C, then by column index.

    Returns (order, scores) where scores is the entry C (inf = never entered).
    Only the first ``n_ranked`` = number of entered columns are a real ranking;
    the remaining tail is in column order and must not be selected.
    """
    X = np.asarray(X, dtype=np.float64)
    std = X.std(axis=0)
    X = (X - X.mean(axis=0)) / np.where(std > 0, std, 1.0)
    classes = np.unique(y)
    entry_step = np.full(X.shape[1], np.inf)
    entry_weight = np.zeros(X.shape[1])
    entry_c = np.full(X.shape[1], np.inf)
    # Slightly below the analytical bound, since class weights shift it.
    C = 0.5 * l1_min_c(X, y, loss="log")
    for index in range(max_steps):
        weights = np.zeros(X.shape[1])
        for label in classes:
            model = LogisticRegression(
                l1_ratio=1.0,
                solver="liblinear",
                C=C,
                class_weight="balanced",
                intercept_scaling=100.0,
                max_iter=5000,
                random_state=seed,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("error", ConvergenceWarning)
                model.fit(X, y == label)
            weights = np.maximum(weights, np.abs(model.coef_[0]))
        new = (weights > 0) & np.isinf(entry_step)
        entry_step[new], entry_weight[new], entry_c[new] = index, weights[new], C
        if np.isfinite(entry_step).sum() >= n_rank:
            break
        C *= step
    else:
        raise RuntimeError(f"Only {np.isfinite(entry_step).sum()} columns entered the L1 path")
    order = np.lexsort((np.arange(X.shape[1]), -entry_weight, entry_step))
    return order, entry_c


RANKERS = {"anova": anova_ranking, "l1": l1_path_ranking}


class RankSelector(SelectorMixin, BaseEstimator):
    """Keep the top ``k`` columns of a training-fold ranking.

    Place it inside a Pipeline so that ``fit`` only ever sees training rows.
    After fitting, ``order_`` holds the full ranking, so any smaller k is a prefix.
    """

    def __init__(self, method="anova", k=50):
        self.method = method
        self.k = k

    def fit(self, X, y):
        X = validate_data(self, X, ensure_all_finite=True)
        if self.method not in RANKERS:
            raise ValueError(f"Unknown selection method: {self.method}")
        if not 1 <= self.k <= X.shape[1]:
            raise ValueError(f"k={self.k} must be between 1 and {X.shape[1]}")
        y = np.asarray(y)
        if self.method == "l1":
            self.order_, self.scores_ = l1_path_ranking(X, y, n_rank=self.k)
        else:
            self.order_, self.scores_ = RANKERS[self.method](X, y)
        # Every method marks unranked columns with an infinite score.
        self.n_ranked_ = int(np.isfinite(self.scores_).sum())
        if self.k > self.n_ranked_:
            raise ValueError(f"Only {self.n_ranked_} columns were ranked; k={self.k}")
        return self

    def _get_support_mask(self):
        check_is_fitted(self, "order_")
        mask = np.zeros(self.n_features_in_, dtype=bool)
        mask[self.order_[: self.k]] = True
        return mask
