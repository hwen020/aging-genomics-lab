"""Shared statistical models and exports for the GSE5086 workbench."""

import json
import platform
from datetime import datetime, timezone
from html import escape
from importlib.metadata import version

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import t

from . import __version__
from .analysis import benjamini_hochberg, qc_summary, sample_pca
from .data import DataValidationError
from .public_studies import AGE_CUTOFF, CARD

MODEL = "Probe-wise OLS with HC3 robust standard errors and two-sided residual-df t tests"
COHORT_COLORS = ["#61D6C4", "#AA9BFF"]
LIMITATIONS = [
    "Exploratory cross-sectional associations. The age-group contrast is defined by age.",
    "Sex is adjusted for. Batch information is unavailable. Other confounding may remain.",
    "The 65-year cutoff is a project choice. The probe-level models do not reproduce the paper's gene-level tests.",
    "Probe mappings can be ambiguous or missing. Results support no causal or clinical inference.",
    "A lack of FDR-significant probes does not establish the absence of an age association.",
]


def fit_probes(ds, predictor, effect_name):
    """Estimate one exposure coefficient across probes with HC3 covariance."""
    metadata = ds.metadata
    sex = metadata.sex.eq("Male").to_numpy(dtype=float)
    design = np.column_stack((np.ones(len(sex)), predictor, sex))
    df = len(sex) - design.shape[1]
    if df <= 0 or np.linalg.matrix_rank(design) != design.shape[1]:
        raise DataValidationError("Rank-deficient sex-adjusted design or insufficient residual degrees of freedom.")
    inverse = np.linalg.pinv(design)
    leverage = np.einsum("ij,ji->i", design, inverse)
    if (leverage >= 1 - 1e-10).any():
        raise DataValidationError("Sample leverage is too high for HC3. Check sex and cohort overlap.")
    outcomes = ds.expression.to_numpy().T
    coefficients = inverse @ outcomes
    residuals = outcomes - design @ coefficients
    se = np.linalg.norm((inverse[1] / (1 - leverage))[:, None] * residuals, axis=0)
    valid = np.isfinite(se) & (se > 1e-12)
    p = np.full(len(se), np.nan)
    fdr = p.copy()
    p[valid] = 2 * t.sf(np.abs(coefficients[1, valid] / se[valid]), df)
    fdr[valid] = benjamini_hochberg(p[valid])
    critical = t.ppf(0.975, df)
    se[~valid] = np.nan
    return pd.DataFrame({
        "gene_id": ds.expression.index, effect_name: coefficients[1],
        "standard_error": se, "ci95_lower": coefficients[1] - critical * se,
        "ci95_upper": coefficients[1] + critical * se,
        "p_value": p, "fdr": fdr,
    })


def analyze_study(study):
    ds, samples, filters, mapping, _ = study
    group = fit_probes(ds, ds.metadata.group.eq("Older"), "log2_fold_change")
    age = fit_probes(ds, (ds.metadata.age - ds.metadata.age.mean()) / 10, "log2_change_per_decade")
    group = group.merge(mapping, on="gene_id", how="left").sort_values(["fdr", "gene_id"])
    age = age.merge(mapping, on="gene_id", how="left").sort_values(["fdr", "gene_id"])
    pca = sample_pca(ds).merge(ds.metadata[["sex"]], left_on="sample_id", right_index=True)
    medians = ds.expression.median(axis=0)
    mad = (medians - medians.median()).abs().median()
    sample_qc = pd.DataFrame({"sample_id": medians.index, "median_log2_signal": medians.values})
    sample_qc["median_robust_z"] = (
        (medians.values - medians.median()) / (1.4826 * mad) if mad > 0 else np.nan
    )
    sample_qc["flag"] = sample_qc.median_robust_z.abs() > 3.5
    return {
        "differential_expression": group, "age_association": age, "sample_pca": pca,
        "sample_selection": samples, "probe_filters": filters, "sample_qc": sample_qc,
        "metadata": ds.metadata.reset_index(),
    }


def make_figures(ds, tables, *, fdr=0.05, effect=0.5, probe=None):
    pca = tables["sample_pca"]
    group = tables["differential_expression"].dropna(subset=["p_value"]).copy()
    age = tables["age_association"].dropna(subset=["p_value"])
    group["Passes thresholds"] = (group.fdr < fdr) & (group.log2_fold_change.abs() >= effect)
    group["−log10 p"] = -np.log10(group.p_value.clip(lower=np.finfo(float).tiny))
    order = ds.metadata.sort_values("age").index
    heat = ds.expression.loc[group.gene_id.head(20), order]
    heat = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1), axis=0)
    charts = {
        "pca": px.scatter(pca, x="pc1", y="pc2", color="group", symbol="sex", hover_data=["sample_id", "age"],
                          color_discrete_sequence=COHORT_COLORS,
                          labels={"pc1": f"PC1 ({pca.pc1_variance_pct.iloc[0]:.1f}%)",
                                  "pc2": f"PC2 ({pca.pc2_variance_pct.iloc[0]:.1f}%)"}),
        "age_distribution": px.box(ds.metadata.reset_index(), x="group", y="age", color="sex", points="all",
                                   hover_data=["sample_id"], color_discrete_sequence=COHORT_COLORS,
                                   labels={"age": "Age (years)", "group": "Age group"}),
        "volcano": px.scatter(group, x="log2_fold_change", y="−log10 p", color="Passes thresholds",
                              hover_data=["gene_id", "gene_symbol", "fdr"],
                              labels={"log2_fold_change": "Sex-adjusted Older − Younger (log2)"}),
        "heatmap": go.Figure(go.Heatmap(z=heat.to_numpy(), x=heat.columns, y=heat.index,
                                        colorscale="RdBu_r", zmid=0, colorbar_title="Row z")),
    }
    charts["heatmap"].update_xaxes(title="GEO samples ordered by age", showticklabels=False)
    if probe is not None or not age.empty:
        probe = probe or age.gene_id.iloc[0]
        points = ds.metadata.assign(log2_signal=ds.expression.loc[probe]).reset_index()
        charts["age_expression"] = px.scatter(
            points, x="age", y="log2_signal", color="sex", hover_data=["sample_id"],
            color_discrete_sequence=COHORT_COLORS,
            title=f"{probe} · observed expression by age",
        )
    return charts


def build_export(study, tables, *, fdr=0.05, effect=0.5, figures=None, age_probe=None):
    """Return the same report files for CLI directories and dashboard ZIPs."""
    ds, samples, filters, _, files = study
    group = tables["differential_expression"]
    age = tables["age_association"]
    qc = qc_summary(ds)
    qc.update({
        "source_samples": len(samples), "source_probes": len(filters),
        "excluded_samples": int((~samples.included).sum()),
        "nonfinite_values_before_probe_filtering": int(filters.nonfinite_values.sum()),
        "zero_values_before_probe_filtering": int(filters.zero_values.sum()),
        "probe_exclusions": filters.loc[filters.exclusion_reason.ne(""), "exclusion_reason"].value_counts().to_dict(),
        "flagged_samples": tables["sample_qc"].loc[tables["sample_qc"].flag, "sample_id"].tolist(),
        "group_tested": int(group.p_value.count()), "age_tested": int(age.p_value.count()),
    })
    provenance = {
        **CARD, "project_version": __version__, "input_type": "public-study", "source_files": files,
        "run_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
        "platform_os": platform.platform(),
        "dependencies": {name: version(name) for name in ("numpy", "pandas", "scipy", "plotly", "streamlit")},
        "random_seed": None, "pca": "Centered probe-wise SVD; no stochastic operations",
        "model": MODEL, "multiple_testing": "BH separately per model across valid tested probes",
        "group_formula": "log2_intensity ~ 1 + I(age >= 65) + I(sex = Male)",
        "age_formula": "log2_intensity ~ 1 + centered_age_years / 10 + I(sex = Male)",
        "covariates": ["sex"], "unavailable_covariates": ["batch"],
        "parameters": {"age_cutoff_years": AGE_CUTOFF, "fdr_strictly_below": fdr,
                       "min_abs_group_log2_effect": effect, "robust_z_flag_above": 3.5,
                       "min_log2_probe_sd_exclusive": 1e-8, "minimum_hc3_se_exclusive": 1e-12,
                       "age_plot_probe": age_probe or next(iter(age.dropna(subset=["p_value"]).gene_id), None)},
        "preprocessing": "Restrict to rectus abdominis; remove nonfinite/zero and near-constant probes; log2 once",
        "qc": qc, "limitations": LIMITATIONS,
    }
    group_hits = int(((group.fdr < fdr) & (group.log2_fold_change.abs() >= effect)).sum())
    age_hits = int((age.fdr < fdr).sum())
    report = f"""# Aging Genomics Lab · {CARD['accession']}

Public human expression microarray · {CARD['citation']}

- GEO record: {CARD['geo_url']}
- Original study: {CARD['paper_url']}
- Access terms: {CARD['access_policy_url']}
- Platform: {CARD['platform']}
- Submitted scale: {CARD['source_scale']}
- Analysis: {qc['n_samples']} rectus abdominis samples from {len(samples)} submitted. {qc['n_features']:,} probe sets retained. See sample and probe logs.
- Scale and preprocessing: {provenance['preprocessing']}
- Model: {MODEL}. Sex-adjusted; batch unavailable.
- Group contrast: Older (≥{AGE_CUTOFF}) minus Younger. {group_hits:,} of {qc['group_tested']:,} tested probes at FDR < {fdr:g} and |log2 effect| ≥ {effect:g}.
- Continuous age: effect per decade. {age_hits:,} of {qc['age_tested']:,} tested probes at FDR < {fdr:g}.
- Multiple testing: {provenance['multiple_testing']}.

## Limitations

""" + "\n".join(f"- {item}" for item in LIMITATIONS) + """

See provenance.json for artifact hashes, acquisition times, formulas, parameters and runtime versions.
"""
    charts = figures if figures is not None else make_figures(ds, tables, fdr=fdr, effect=effect)
    artifacts = {f"{name}.csv": table.to_csv(index=False) for name, table in tables.items()}
    artifacts["report.md"] = report
    artifacts["report.html"] = '<!doctype html><meta charset="utf-8"><title>GSE5086 report</title><pre>' + escape(report) + "</pre>"
    artifacts["provenance.json"] = json.dumps(provenance, indent=2) + "\n"
    for name, chart in charts.items():
        artifacts[f"chart_{name}.html"] = chart.to_html(include_plotlyjs=True)
    return artifacts, provenance
