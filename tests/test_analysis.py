import numpy as np
import pytest

from aginggenomics.analysis import (age_association, benjamini_hochberg,
                                    differential_expression, qc_summary, sample_pca)
from aginggenomics.data import validate_dataset
from aginggenomics.simulation import make_demo


@pytest.fixture(scope="module")
def ds():
    expr, meta, _ = make_demo()
    return validate_dataset(expr, meta, simulated=True)


def test_bh_example():
    got = benjamini_hochberg([0.01, 0.04, 0.03, 0.002])
    np.testing.assert_allclose(got, [0.02, 0.04, 0.04, 0.008])


def test_bh_preserves_monotonicity_and_bounds():
    p = np.array([.8, .01, .1, .03, .04, 1, 0])
    q = benjamini_hochberg(p)
    assert np.all((q >= 0) & (q <= 1))
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)


def test_bh_bad_p_rejected():
    with pytest.raises(ValueError):
        benjamini_hochberg([-0.2, 1.1])


def test_bh_empty():
    assert len(benjamini_hochberg([])) == 0


def test_de_recovers_injected_signal(ds):
    de = differential_expression(ds, "Younger", "Older")
    assert len(de) == 320
    assert de["gene_id"].is_unique
    assert (de["fdr"] >= 0).all()
    assert de.set_index("gene_id").loc["SIM_GENE_0001", "log2_fold_change"] > .7
    assert de.set_index("gene_id").loc["SIM_GENE_0014", "log2_fold_change"] < -.7
    injected = [f"SIM_GENE_{i:04d}" for i in range(1, 25)]
    assert (de.set_index("gene_id").loc[injected, "fdr"] < 0.05).sum() >= 17


def test_de_reverse_contrast_flips_effect(ds):
    a = differential_expression(ds, "Younger", "Older").set_index("gene_id")
    b = differential_expression(ds, "Older", "Younger").set_index("gene_id")
    np.testing.assert_allclose(a.log2_fold_change, -b.loc[a.index].log2_fold_change)


def test_de_invalid_contrast(ds):
    with pytest.raises(ValueError, match="distinct"):
        differential_expression(ds, "Younger", "Younger")


def test_age_association_has_valid_p(ds):
    age = age_association(ds)
    assert len(age) == 320
    assert age["fdr"].between(0, 1).all()
    assert age["spearman_rho"].between(-1, 1).all()


def test_age_assoc_subgroup(ds):
    age = age_association(ds, group="Younger")
    assert age.shape == (320, 4)


def test_age_assoc_unknown_group(ds):
    with pytest.raises(ValueError, match="Unknown"):
        age_association(ds, group="Missing")


def test_pca_shape_and_variance(ds):
    pca = sample_pca(ds)
    assert len(pca) == ds.expression.shape[1]
    assert pca.pc1_variance_pct.iloc[0] >= pca.pc2_variance_pct.iloc[0]
    assert 0 < pca.pc1_variance_pct.iloc[0] + pca.pc2_variance_pct.iloc[0] <= 100.00001


def test_qc(ds):
    qc = qc_summary(ds)
    assert qc["n_samples"] == 60
    assert qc["n_features"] == 320
    assert qc["nonfinite_values"] == 0
    assert abs(qc["age_group_correlation"]) > .8
