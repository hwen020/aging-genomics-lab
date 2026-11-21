"""Independent numerical checks for the two public statistical models."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.special import betainc, expit

from aginggenomics.analysis import benjamini_hochberg
from aginggenomics.data import Dataset, DataValidationError
from aginggenomics.methylation_analysis import age_design, fit_block
from aginggenomics.public_analysis import fit_probes


def reference_hc3(x, y):
    """Scalar sandwich covariance and Student tail via the incomplete beta function."""
    bread = np.linalg.inv(x.T @ x)
    b = np.linalg.solve(x.T @ x, x.T @ y)
    meat = np.zeros_like(bread)
    for row, observed in zip(x, y):
        residual = observed - row @ b
        h = row @ bread @ row
        meat += np.outer(row, row) * (residual / (1 - h)) ** 2
    se = np.sqrt((bread @ meat @ bread)[1, 1])
    df = len(y) - x.shape[1]
    p = betainc(df / 2, 0.5, df / (df + (b[1] / se) ** 2))
    return b[1], se, p


def test_bh_known_values_with_ties_and_unsorted_order():
    np.testing.assert_allclose(benjamini_hochberg([0.04, 0, 0.01, 0.04, 1]),
                               [0.05, 0, 0.025, 0.05, 1])


@pytest.mark.parametrize("group_effect", [False, True])
def test_expression_hc3_signed_effect_matches_independent_reference(group_effect):
    age = np.array([25, 32, 47, 55, 67, 71, 80, 85])
    male = np.array([0, 1] * 4)
    predictor = (age >= 65) if group_effect else (age - age.mean()) / 10
    design = np.column_stack([np.ones(8), predictor, male])
    y = 7 - 0.37 * predictor + 0.2 * male + np.array([.02, -.04, .01, .06, -.05, .04, .02, -.03])
    metadata = pd.DataFrame({"age": age, "sex": np.where(male, "Male", "Female")})
    ds = Dataset(pd.DataFrame([y, -y], index=["probe-1", "probe-2"]), metadata)
    result = fit_probes(ds, predictor, "known_effect")
    effect, se, p = reference_hc3(design, y)
    assert effect < 0
    np.testing.assert_allclose(result.known_effect, [effect, -effect], rtol=1e-11)
    np.testing.assert_allclose(result.standard_error, [se, se], rtol=1e-11)
    np.testing.assert_allclose(result.p_value, [p, p], rtol=1e-10)


def test_ewas_hc3_boundary_clipping_and_plate_reference():
    age = np.array([19, 28, 37, 43, 59, 64, 70, 88, 30, 51, 74, 95])
    male = np.array([0, 1] * 6)
    plate = np.array(["1", "10", "2", "1"] * 3)
    x = np.column_stack([np.ones(len(age)), (age - age.mean()) / 10, male, plate == "10", plate == "2"])
    meta = pd.DataFrame({"age": age, "sex": np.where(male, "M", "F"), "plate": plate})
    model = age_design(meta)
    np.testing.assert_array_equal(model[0], x)
    beta = expit(np.log(2) * (-0.2 + 0.65 * x[:, 1] - 0.4 * male + x[:, 3] + 0.3 * x[:, 4]))
    beta[0], beta[-1] = 0, 1
    bounded = [min(1 - 1e-6, max(1e-6, b)) for b in beta]
    y = np.array([(np.log(b) - np.log1p(-b)) / np.log(2) for b in bounded])
    actual = np.array(fit_block(beta[None, :], model)).ravel()
    expected = reference_hc3(x, y)
    assert expected[0] > 0
    np.testing.assert_allclose(actual, expected, rtol=1e-10)
    np.testing.assert_allclose(fit_block((1 - beta)[None, :], model)[0], [-expected[0]], rtol=1e-10)
    with pytest.raises(DataValidationError, match="Rank-deficient"):
        age_design(meta.assign(plate=np.where(male, "1", "2")))


@pytest.mark.parametrize(("accession", "tests", "hits"), [("GSE5086", 54615, 0), ("GSE40279", 424870, 58849)])
def test_portfolio_is_read_only_and_reports_verified_counts(monkeypatch, accession, tests, hits):
    from streamlit.testing.v1 import AppTest

    monkeypatch.setenv("AGING_GENOMICS_MODE", "portfolio")
    app = AppTest.from_file(Path(__file__).resolve().parents[1] / "app.py", default_timeout=60).run()
    assert not app.get("file_uploader")
    app.sidebar.radio[0].set_value("Public Study").run()
    app.selectbox[0].set_value(accession).run()
    assert not app.exception and not app.error
    assert app.info[0].value == "PUBLIC STUDY · PRECOMPUTED SNAPSHOT"
    assert [item.value for item in app.metric][1:] == [f"{tests:,}", f"{hits:,}"]
    assert not app.button and len(app.get("plotly_chart")) == 1
