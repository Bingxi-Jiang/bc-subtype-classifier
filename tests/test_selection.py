import numpy as np
import pytest
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline

from bc_subtype.selection import RankSelector, anova_ranking, l1_path_ranking


def toy_data():
    y = np.array(["a"] * 4 + ["b"] * 4)
    shift = np.repeat([0.0, 1.0], 4)
    noise = np.array([0.1, -0.1, 0.2, -0.2, 0.1, -0.1, 0.2, -0.2])
    X = np.column_stack(
        [
            noise,  # 0: no class signal
            shift * 5 + noise,  # 1: strongest signal
            shift * 1 + noise,  # 2: weak signal
            np.full(8, 3.0),  # 3: constant -> undefined F
            shift * 5 + noise,  # 4: exact tie with column 1
        ]
    )
    return X, y


def test_anova_ranking_orders_by_signal_breaks_ties_by_index_and_puts_constants_last():
    X, y = toy_data()
    order, scores = anova_ranking(X, y)
    assert order.tolist() == [1, 4, 2, 0, 3]
    assert scores[3] == -np.inf
    assert sorted(order.tolist()) == list(range(X.shape[1]))


def test_selector_keeps_nested_top_k_and_validates_inputs():
    X, y = toy_data()
    masks = {k: RankSelector("anova", k).fit(X, y).get_support() for k in (1, 2, 3)}
    assert np.flatnonzero(masks[2]).tolist() == [1, 4]
    assert all((masks[k] <= masks[k + 1]).all() for k in (1, 2))
    assert RankSelector("anova", 2).fit(X, y).transform(X).shape == (8, 2)
    with pytest.raises(ValueError, match="must be between"):
        RankSelector("anova", 6).fit(X, y)
    with pytest.raises(ValueError, match="Unknown selection"):
        RankSelector("nope", 2).fit(X, y)


def test_selector_in_pipeline_sees_only_training_rows():
    X, y = toy_data()
    train = np.array([0, 1, 2, 4, 5, 6])
    # Make column 0 informative only in rows outside the training set; a selector
    # that peeked at those rows would rank it differently.
    X_leak = X.copy()
    X_leak[[3, 7], 0] = [-100.0, 100.0]
    model = Pipeline([("select", RankSelector("anova", 2)), ("clf", LogisticRegression())])
    model.fit(X_leak[train], y[train])
    reference = RankSelector("anova", 2).fit(X[train], y[train])
    np.testing.assert_array_equal(model.named_steps["select"].order_, reference.order_)


def redundant_data():
    # Columns 0 and 1 carry the same class-"a" signal; column 2 is the only one that
    # separates "b" from "c"; columns 3+ are noise.
    rng = np.random.default_rng(0)
    y = np.repeat(["a", "b", "c"], 30)
    signal_a = (y == "a") * 3.0 + rng.normal(0, 0.5, 90)
    X = np.column_stack(
        [
            signal_a,
            signal_a + rng.normal(0, 0.05, 90),
            (y == "c") * 1.5 + rng.normal(0, 0.5, 90),
            rng.normal(0, 1, (90, 7)),
        ]
    )
    return X, y


def test_l1_path_prefers_complementary_over_redundant_columns():
    X, y = redundant_data()
    anova_order, _ = anova_ranking(X, y)
    assert set(anova_order[:2]) == {0, 1}  # univariate: both copies first
    order, entry_c = l1_path_ranking(X, y, n_rank=3)
    assert 2 in order[:2] and len({0, 1} & set(order[:2])) == 1
    ranked = np.isfinite(entry_c)
    assert ranked.sum() >= 3
    assert np.all(np.diff(entry_c[order[: ranked.sum()]]) >= 0)  # entry C never decreases
    assert np.isinf(entry_c[order[ranked.sum() :]]).all()


def test_l1_selector_and_path_limits():
    X, y = redundant_data()
    selector = RankSelector("l1", k=3).fit(X, y)
    assert selector.get_support().sum() == 3 and selector.n_ranked_ >= 3
    with pytest.raises(RuntimeError, match="entered the L1 path"):
        l1_path_ranking(X, y, n_rank=10, max_steps=1)
