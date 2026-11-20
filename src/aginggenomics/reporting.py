"""Report and provenance artifacts for reproducible local analysis."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from .analysis import age_association, differential_expression, qc_summary, sample_pca
from .data import load_dataset
from . import __version__


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def export_analysis(expression_path: str | Path, metadata_path: str | Path,
                    out_dir: str | Path, *, reference: str | None = None,
                    target: str | None = None, simulated: bool = False,
                    source: str = "user-imported") -> dict:
    ds = load_dataset(expression_path, metadata_path, source=source, simulated=simulated)
    ref = reference or ("Younger" if "Younger" in ds.groups else ds.groups[0])
    tgt = target or next(g for g in ds.groups if g != ref)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    de = differential_expression(ds, ref, tgt)
    age = age_association(ds)
    pca = sample_pca(ds)
    qc = qc_summary(ds)
    de.to_csv(out / "differential_expression.csv", index=False)
    age.to_csv(out / "age_association.csv", index=False)
    pca.to_csv(out / "sample_pca.csv", index=False)
    provenance = {
        "project_version": __version__,
        "analysis_time_utc": datetime.now(timezone.utc).isoformat(),
        "input_type": "synthetic" if simulated else "user-imported",
        "source": source,
        "expression_sha256": sha256(expression_path),
        "metadata_sha256": sha256(metadata_path),
        "reference_group": ref, "target_group": tgt,
        "statistical_method": "Welch t-test on normalized log2 expression + BH FDR",
        "age_method": "Unadjusted Spearman correlation + BH FDR",
        "qc": qc,
        "limitations": [
            "Not suitable for raw RNA-seq counts",
            "Age and groups may be confounded; no covariate adjustment",
            "Not for clinical, diagnostic, or causal claims",
        ],
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    hits_de = de[(de.fdr < 0.05) & (de.log2_fold_change.abs() >= 0.5)]
    hits_age = age[age.fdr < 0.05]
    banner = "SYNTHETIC DEMONSTRATION — NOT REAL BIOLOGICAL FINDINGS" if simulated else "EXPLORATORY ANALYSIS — VERIFY PREPROCESSING AND STUDY DESIGN"
    report = f"""# Aging Genomics — Analysis Report

> **{banner}**

- Source: `{source}`
- Data: {qc['n_samples']} samples; {qc['n_features']} features
- Contrast: {tgt} vs {ref}
- Expression interpretation: **pre-normalized log2 values**
- Differential method: Welch t-tests, BH-FDR; |log2 difference| ≥ 0.5
- Age association: Spearman correlation, BH-FDR; **unadjusted for covariates**
- Significant group features (FDR < 0.05 and effect ≥ 0.5): **{len(hits_de)}**
- Age-associated features (FDR < 0.05): **{len(hits_age)}**

## Interpretation limitations

These are exploratory statistics; age-group confounding, batch effects, tissue,
sex and other covariates can invalidate biological interpretation.
Raw RNA-seq counts need a count-aware model, not this Welch t-test workflow.
No result here is a medical or clinical inference.

## Reproducibility

`provenance.json` contains input SHA-256 hashes, methods, parameters and the
software version. Tables can be regenerated with the documented CLI command.
"""
    (out / "report.md").write_text(report, encoding="utf-8")
    return provenance
