"""Standalone, reproducible EWAS reports and interactive Plotly figures."""

import json
import platform
import tempfile
from datetime import datetime, timezone
from html import escape
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px

from . import __version__
from .methylation_analysis import CAVEATS, FDR
from .methylation_study import ACCESSION, GEO, PAPER, TERMS


def figures(results, pca, betas, *, probe=None):
    tested = results.loc[results.p_value.notna()]
    ordered = tested.sort_values("cpg_id")
    step = max(1, (len(ordered) + 8999) // 9000)
    volcano = pd.concat([tested.head(1000), ordered.iloc[::step]]).drop_duplicates("cpg_id")
    volcano = volcano.assign(significant=volcano.fdr.lt(FDR),
                             neg_log10_p=-np.log10(volcano.p_value.clip(lower=np.finfo(float).tiny)))
    probe = probe or next(iter(betas.cpg_id.unique()))
    trace = betas.loc[betas.cpg_id.eq(probe)]
    pca_figure = px.scatter(pca, x="pc1", y="pc2", color="age", symbol="sex",
                           hover_data=["sample_id", "plate"],
                           labels={"age": "Age (years)", "pc1": f"PC1 ({pca.pc1_variance_pct.iloc[0]:.1f}%)",
                                   "pc2": f"PC2 ({pca.pc2_variance_pct.iloc[0]:.1f}%)"})
    pca_figure.update_layout(legend={"orientation": "h", "x": 0, "y": 1.04})
    return {
        "pca": pca_figure,
        "volcano": px.scatter(volcano, x="m_change_per_decade", y="neg_log10_p", color="significant",
                              hover_data=["cpg_id", "gene_symbols", "fdr"],
                              labels={"m_change_per_decade": "Adjusted M-value change per decade",
                                      "neg_log10_p": "−log10 p", "significant": "FDR < 0.05"}),
        "age_beta": px.scatter(trace, x="age", y="beta", color="sex",
                               hover_data=["sample_id", "plate"], title=f"{probe} · submitted beta by age",
                               labels={"age": "Age (years)", "beta": "Observed beta ([0, 1])"}),
    }


def export_ewas(run, inputs, out):
    results = run.results
    hits = results.loc[results.fdr < FDR]
    filters = run.probe_qc.exclusion_reason
    qc = {
        "samples": len(run.metadata), "age_range": [int(run.metadata.age.min()), int(run.metadata.age.max())],
        "source_probes": len(filters), "eligible_probes": len(results),
        "tested_probes": int(results.p_value.count()), "fdr_hits": len(hits),
        "sample_outlier_flags": int(run.sample_qc.flag_outlier.sum()),
        "sample_missingness_flags": int(run.sample_qc.flag_missing.sum()),
        "excluded_probes": filters[filters.ne("")].value_counts().to_dict(),
    }
    provenance = {
        "accession": ACCESSION, "geo_url": GEO, "paper_url": PAPER, "terms_url": TERMS,
        "downloaded_artifacts": inputs, "run_utc": datetime.now(timezone.utc).isoformat(),
        "project_version": __version__, "python": platform.python_version(), "system": platform.platform(),
        "dependencies": {name: version(name) for name in ("numpy", "pandas", "scipy", "plotly", "streamlit")},
        "data_origin": "Public Study", "platform": "GPL13534 Illumina HumanMethylation450 hg19",
        "submitted_scale": "Author-processed beta fractions [0, 1]", "model_scale": "log2 beta odds",
        "inference": "OLS with HC3 robust standard errors, two-sided t tests, 95% intervals and BH",
        "covariates": ["sex", "plate"], "sex_reference": "F",
        "plate_reference": sorted(run.metadata.plate.unique())[0],
        "parameters": run.parameters, "qc": qc, "limitations": CAVEATS,
    }
    report = f"""# GSE40279 · Continuous-age EWAS

Public whole-blood methylation study on GPL13534. Hannum et al., Molecular Cell (2013).
[GEO]({GEO}) · [Paper]({PAPER}) · [Access terms]({TERMS})

Analyzed **{qc['samples']} samples** aged {qc['age_range'][0]}–{qc['age_range'][1]}.
Of {qc['source_probes']:,} input probes, {qc['tested_probes']:,} had valid age tests
and **{qc['fdr_hits']:,} CpGs** passed FDR < {FDR}. See probe and sample QC tables.

## Methods

Submitted beta fractions were checked against [0, 1]. Tests use M-values
computed as log2 odds after clipping beta at {run.parameters['beta_clip_epsilon']:g}
and its upper complement. Probe-wise OLS estimates change per decade adjusted for
recorded sex and categorical plate. HC3 robust standard errors give two-sided
t tests with {run.parameters['residual_df']} residual degrees of freedom and
95% confidence intervals. BH covers all valid tests in this run.

Probe QC excludes missing/out-of-range values, non-CpG and non-autosomal probes,
incomplete hg19 annotations, probes with a GPL13534 SNP annotation within 10 bases
and beta SD below {run.parameters['min_beta_sd']:g}. Sample mean-beta and missingness
flags are descriptive and all matched samples are retained.
Descriptive PCA uses {run.parameters['pca_probes']:,} eligible beta probes selected
deterministically by CpG ID modulo {run.parameters['pca_probe_id_modulus']}.
Three interactive chart files accompany the full table. The volcano shows the top
1,000 tests plus a deterministic subset of up to 9,000 tested probes.

## Limitations

""" + "\n".join(f"- {item}" for item in CAVEATS) + "\n\nSee provenance.json for checksums, design, parameters and versions.\n"
    betas = pd.concat([run.metadata[["sample_id", "age", "sex", "plate"]].assign(cpg_id=cpg, beta=beta)
                       for cpg, beta in run.example_betas.items()], ignore_index=True)
    charts = figures(results, run.pca, betas)
    tables = {"ewas_results": results, "significant_sites": hits, "sample_metadata": run.metadata,
              "probe_qc": run.probe_qc, "sample_qc": run.sample_qc, "sample_pca": run.pca, "top_cpg_beta": betas}
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out) as staging:
        folder = Path(staging)
        for name, table in tables.items():
            table.to_csv(folder / f"{name}.csv", index=False)
        for name, figure in charts.items():
            figure.write_html(folder / f"chart_{name}.html", include_plotlyjs=True)
        (folder / "report.md").write_text(report, encoding="utf-8")
        (folder / "report.html").write_text(
            '<!doctype html><meta charset="utf-8"><title>GSE40279 EWAS report</title>'
            '<main style="max-width:78ch;margin:3rem auto;font:16px system-ui">'
            '<pre style="white-space:pre-wrap;font:inherit">' + escape(report) + "</pre></main>", encoding="utf-8")
        (folder / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
        for item in sorted(folder.iterdir(), key=lambda p: p.name == "provenance.json"):
            item.replace(out / item.name)
    return provenance
