"""Transparent, non-clinical exploratory statistics for normalized log2 expression."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_ind

from .data import Dataset


def benjamini_hochberg(pvalues: np.ndarray | list[float]) -> np.ndarray:
    """Benjamini-Hochberg false-discovery-rate adjusted p-values."""
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be a finite vector in [0, 1].")
    if len(p) == 0:
        return np.array([], dtype=float)
    ordering = np.argsort(p)
    ranked = p[ordering]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[ordering] = np.minimum(adjusted, 1.0)
    return result


def differential_expression(ds: Dataset, reference: str = "Younger",
                            target: str = "Older") -> pd.DataFrame:
    """Per-gene Welch tests. Valid only on properly pre-normalized log2 expression.

    log2FC = mean(target) - mean(reference). Not a raw-count RNA-seq DE method.
    """
    if {reference, target} != set(ds.groups) or reference == target:
        raise ValueError(f"Select two distinct groups from {ds.groups}.")
    refs = ds.metadata.index[ds.metadata.group == reference]
    tgts = ds.metadata.index[ds.metadata.group == target]
    a = ds.expression[refs].to_numpy(dtype=float)
    b = ds.expression[tgts].to_numpy(dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        pvalue = ttest_ind(b, a, equal_var=False, axis=1).pvalue
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
    pvalue = np.clip(pvalue, 0.0, 1.0)
    mean_ref = np.mean(a, axis=1)
    mean_target = np.mean(b, axis=1)
    result = pd.DataFrame({
        "gene_id": ds.expression.index.to_numpy(),
        "mean_reference": mean_ref,
        "mean_target": mean_target,
        "log2_fold_change": mean_target - mean_ref,
        "p_value": pvalue,
        "fdr": benjamini_hochberg(pvalue),
    })
    return result.sort_values(["fdr", "p_value", "gene_id"], kind="stable").reset_index(drop=True)


def age_association(ds: Dataset, *, group: str | None = None) -> pd.DataFrame:
    """Exploratory per-gene Spearman associations with chronological age.

    P-values do not adjust for sex, batch, ancestry, cell composition or other
    covariates. Group-vs-age confounding is especially important in aging cohorts.
    """
    if group is not None and group not in ds.groups:
        raise ValueError(f"Unknown group: {group}")
    samples = ds.metadata.index if group is None else ds.metadata.index[ds.metadata.group == group]
    ages = ds.metadata.loc[samples, "age"].to_numpy(dtype=float)
    if len(np.unique(ages)) < 3:
        raise ValueError("At least three distinct ages are needed for association analysis.")
    rows = []
    for gene, vals in ds.expression[samples].iterrows():
        x = vals.to_numpy(dtype=float)
        if np.std(x) <= 1e-12:
            rho, p = 0.0, 1.0
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                rho, p = spearmanr(ages, x)
            if not np.isfinite(p):
                rho, p = 0.0, 1.0
        rows.append((gene, float(rho), float(p)))
    frame = pd.DataFrame(rows, columns=["gene_id", "spearman_rho", "p_value"])
    frame["fdr"] = benjamini_hochberg(frame["p_value"].to_numpy())
    return frame.sort_values(["fdr", "p_value", "gene_id"], kind="stable").reset_index(drop=True)


def sample_pca(ds: Dataset) -> pd.DataFrame:
    """Centered-feature PCA by SVD; no test-label leakage or fitting to target."""
    matrix = ds.expression.T.to_numpy(dtype=float)
    matrix = matrix - matrix.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    coords = u[:, :2] * s[:2]
    weights = s ** 2
    var_ratio = weights[:2] / weights.sum() if weights.sum() > 0 else np.zeros(2)
    return pd.DataFrame({
        "sample_id": ds.metadata.index.to_list(),
        "pc1": coords[:, 0],
        "pc2": coords[:, 1],
        "group": ds.metadata["group"].to_list(),
        "age": ds.metadata["age"].to_list(),
        "pc1_variance_pct": float(var_ratio[0] * 100),
        "pc2_variance_pct": float(var_ratio[1] * 100),
    })


def qc_summary(ds: Dataset) -> dict:
    values = ds.expression.to_numpy()
    return {
        "n_samples": ds.expression.shape[1],
        "n_features": ds.expression.shape[0],
        "group_counts": ds.metadata.group.value_counts().to_dict(),
        "age_min": float(ds.metadata.age.min()),
        "age_max": float(ds.metadata.age.max()),
        "nonfinite_values": int((~np.isfinite(values)).sum()),
        "near_constant_features": int((np.std(values, axis=1) < 1e-8).sum()),
        "age_group_correlation": float(np.corrcoef(
            ds.metadata["age"].to_numpy(),
            (ds.metadata["group"] == ds.groups[1]).astype(float).to_numpy()
        )[0, 1]),
    }
