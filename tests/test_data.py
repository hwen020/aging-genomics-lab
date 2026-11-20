from io import StringIO

import numpy as np
import pandas as pd
import pytest

from aginggenomics.data import DataValidationError, load_dataset, validate_dataset
from aginggenomics.simulation import make_demo


@pytest.fixture
def demo():
    return make_demo(n_genes=100, per_group=8)


def test_seed_reproducible():
    x, y, z = make_demo()
    a, b, c = make_demo()
    pd.testing.assert_frame_equal(x, a)
    pd.testing.assert_frame_equal(y, b)
    pd.testing.assert_frame_equal(z, c)


def test_group_aligns_samples_after_metadata_shuffle(demo):
    expr, meta, _ = demo
    ds = validate_dataset(expr, meta.sample(frac=1, random_state=1), simulated=True)
    assert ds.expression.columns.tolist() == ds.metadata.index.tolist()
    assert ds.groups == ["Older", "Younger"]


def test_missing_metadata_rejected(demo):
    expr, meta, _ = demo
    with pytest.raises(DataValidationError, match="Sample IDs differ"):
        validate_dataset(expr, meta.iloc[:-1])


def test_duplicate_gene_rejected(demo):
    expr, meta, _ = demo
    expr.loc[1, "gene_id"] = expr.loc[0, "gene_id"]
    with pytest.raises(DataValidationError, match="unique"):
        validate_dataset(expr, meta)


def test_non_numeric_rejected(demo):
    expr, meta, _ = demo
    expr[expr.columns[1]] = expr[expr.columns[1]].astype(object)
    expr.iloc[0, 1] = "invalid"
    with pytest.raises(DataValidationError, match="nonnumeric"):
        validate_dataset(expr, meta)


def test_missing_value_rejected(demo):
    expr, meta, _ = demo
    expr.iloc[0, 1] = np.nan
    with pytest.raises(DataValidationError, match="missing"):
        validate_dataset(expr, meta)


def test_raw_counts_obviously_too_large_rejected(demo):
    expr, meta, _ = demo
    expr.iloc[:, 1:] = 1000
    with pytest.raises(DataValidationError, match="raw counts"):
        validate_dataset(expr, meta)


def test_small_integer_count_like_values_rejected(demo):
    expr, meta, _ = demo
    expr.iloc[:, 1:] = expr.iloc[:, 1:].round().abs()
    with pytest.raises(DataValidationError, match="ambiguous with raw counts"):
        validate_dataset(expr, meta)


def test_group_small_rejected(demo):
    expr, meta, _ = demo
    meta.loc[meta.group == "Older", "group"] = "Younger"
    with pytest.raises(DataValidationError, match="Exactly two"):
        validate_dataset(expr, meta)


def test_csv_files_read(demo):
    expr, meta, _ = demo
    ds = load_dataset(StringIO(expr.to_csv(index=False)), StringIO(meta.to_csv(index=False)))
    assert ds.expression.shape == (100, 16)
    assert not ds.is_simulated


def test_invalid_age_rejected(demo):
    expr, meta, _ = demo
    meta.loc[0, "age"] = 999
    with pytest.raises(DataValidationError, match="0–120"):
        validate_dataset(expr, meta)


def test_insufficient_unique_ages_rejected(demo):
    expr, meta, _ = demo
    meta["age"] = 50
    with pytest.raises(DataValidationError, match="distinct ages"):
        validate_dataset(expr, meta)
