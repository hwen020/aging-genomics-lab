"""Aging Genomics Lab — local Streamlit research-software portfolio demo."""
from __future__ import annotations

from io import BytesIO
import hashlib
import json
import os
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from aginggenomics import __version__
from aginggenomics.analysis import age_association, differential_expression, qc_summary, sample_pca
from aginggenomics.data import DataValidationError, load_dataset


st.set_page_config(page_title="Aging Genomics Lab", page_icon="🧬", layout="wide",
                   initial_sidebar_state="expanded")
st.markdown("""
<style>
.block-container {max-width: 1480px; padding-top: 2rem; padding-bottom: 4rem;}
h1, h2, h3 {letter-spacing: -.035em;}
[data-testid="stSidebar"] {border-right: 1px solid #283348;}
[data-testid="stMetric"] {background:#141e31; border:1px solid #253348;
  border-radius:14px; padding:17px 18px;}
[data-testid="stMetricLabel"] {color:#A5B4C8;}
[data-testid="stMetricValue"] {font-weight:700;}
.hero {padding:28px 30px 24px;border:1px solid #29374c;border-radius:20px;
 background:radial-gradient(circle at 96% 4%,#164a55 0,transparent 42%),
 linear-gradient(100deg,#132237,#0c1424);margin-bottom:24px;}
.eyebrow {color:#74e6ce;letter-spacing:.15em;font-size:.72rem;font-weight:700;}
.hero-title {color:#F2F7FC;font-size:2.35rem;font-weight:740;letter-spacing:-.05em;margin:10px 0 7px;}
.hero-sub {color:#B2C2D5;font-size:.99rem;max-width:850px;line-height:1.62;}
.muted {color:#95a9c1;font-size:.83rem;}
.warning-panel {padding:12px 16px;border:1px solid #80653c;background:#332b20;
 border-radius:10px;color:#f8ddb0;margin:6px 0 15px;line-height:1.6;}
.kicker {font-size:.78rem;color:#8ea5c0;text-transform:uppercase;letter-spacing:.09em;}
</style>""", unsafe_allow_html=True)

PLOT_BG = "rgba(0,0,0,0)"
GRID = "#283448"
TEAL = "#61D6C4"
PURPLE = "#AA9BFF"


def style_fig(fig: go.Figure, height: int = 420) -> go.Figure:
    fig.update_layout(template="plotly_dark", height=height, paper_bgcolor=PLOT_BG,
                      plot_bgcolor=PLOT_BG, font=dict(color="#C6D3E4", size=12),
                      margin=dict(t=32, b=35, l=40, r=22), legend_title_text="")
    fig.update_xaxes(gridcolor=GRID, zerolinecolor=GRID)
    fig.update_yaxes(gridcolor=GRID, zerolinecolor=GRID)
    return fig


def csv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False).encode("utf-8")


def archive_bytes(ds, de, age, pca, qc, ref, tgt, input_hashes) -> bytes:
    bio = BytesIO()
    manifest = {
        "project_version": __version__,
        "data_origin": "synthetic" if ds.is_simulated else "user-imported",
        "source": ds.source,
        "reference_group": ref,
        "target_group": tgt,
        "input_sha256": input_hashes,
        "qc": qc,
        "statistics": "Welch t-test on pre-normalized log2 expression; BH FDR",
        "age_association": "Unadjusted Spearman; BH FDR",
        "limitations": "Exploratory only; raw counts unsupported; confounding not adjusted",
    }
    report = ("# Aging Genomics Lab — Export\n\n"
              f"Dataset: **{'SIMULATED, not human' if ds.is_simulated else 'User-imported'}**\n\n"
              f"Contrast: {tgt} vs {ref}; samples: {qc['n_samples']}; "
              f"features: {qc['n_features']}\n\n"
              "**Methods:** Welch t-test on pre-normalized log2 expression, "
              "Benjamini–Hochberg FDR. Spearman age association is unadjusted.\n\n"
              "**Caution:** No clinical interpretation; do not analyze raw RNA-seq "
              "counts with this model. Group and age may be confounded.\n")
    with ZipFile(bio, "w", ZIP_DEFLATED) as z:
        z.writestr("results/differential_expression.csv", csv_bytes(de))
        z.writestr("results/age_association.csv", csv_bytes(age))
        z.writestr("results/sample_pca.csv", csv_bytes(pca))
        z.writestr("results/report.md", report)
        z.writestr("results/provenance.json", json.dumps(manifest, indent=2))
    return bio.getvalue()


portfolio = os.getenv("AGING_GENOMICS_MODE", "research") != "research"
if portfolio:
    st.sidebar.caption("Portfolio Demo · public snapshots")
mode = st.sidebar.radio("Workspace", ["Synthetic Demo" if portfolio else "Synthetic Demo / CSV upload", "Public Study"])
if mode == "Public Study":
    if portfolio:
        from aginggenomics.portfolio_dashboard import render_portfolio

        render_portfolio(style_fig)
    else:
        from aginggenomics.public_dashboard import render_public_study

        render_public_study(style_fig)
    st.stop()


with st.sidebar:
    st.markdown("### 🧬 AGING GENOMICS")
    st.caption("RESEARCH SOFTWARE · v" + __version__)
    st.divider()
    uploaded_expr = uploaded_meta = None
    if not portfolio:
        st.subheader("Dataset")
        st.caption("Upload **both** CSV files for your own pre-normalized log2 expression data.")
        uploaded_expr = st.file_uploader("Expression matrix", type="csv", key="expression")
        uploaded_meta = st.file_uploader("Sample metadata", type="csv", key="metadata")
        st.caption("Expected: gene_id × sample columns, and sample_id / group / age metadata.")
    st.divider()
    st.markdown("**Analysis defaults**")
    fdr_threshold = st.slider("FDR cutoff", 0.01, 0.20, 0.05, 0.01)
    effect_threshold = st.slider("Minimum |log2 difference|", 0.0, 2.0, 0.5, 0.1)
    st.divider()
    if not portfolio:
        st.caption("No user uploads are saved by this app. Local hosting only unless you choose to deploy it.")

if bool(uploaded_expr) != bool(uploaded_meta):
    st.info("Upload the expression matrix **and** metadata file to switch from the simulated demo.")

try:
    if uploaded_expr is not None and uploaded_meta is not None:
        with st.sidebar:
            input_type = st.selectbox("Expression data type", ["Normalized log2 expression", "Raw RNA-seq counts"])
            if input_type == "Raw RNA-seq counts":
                st.error("Raw RNA-seq counts require a count-aware model such as DESeq2.")
                st.stop()
            confirmed_scale = st.checkbox(
                "I confirm that values are pre-normalized log2 expression (not raw counts).",
                value=False)
        if not confirmed_scale:
            st.warning("Confirm the input expression scale in the sidebar before analyzing imported data.")
            st.stop()
        raw_expr = uploaded_expr.getvalue()
        raw_meta = uploaded_meta.getvalue()
        ds = load_dataset(BytesIO(raw_expr), BytesIO(raw_meta),
                          source=f"uploaded: {uploaded_expr.name} + {uploaded_meta.name}")
    else:
        with open("data/demo/expression.csv", "rb") as f:
            raw_expr = f.read()
        with open("data/demo/metadata.csv", "rb") as f:
            raw_meta = f.read()
        ds = load_dataset(BytesIO(raw_expr), BytesIO(raw_meta),
                          source="SYNTHETIC_AGING_DEMO_V1", simulated=True)
except (DataValidationError, OSError) as exc:
    st.error(f"Dataset validation failed: {exc}")
    st.stop()

input_hashes = {"expression": hashlib.sha256(raw_expr).hexdigest(),
                "metadata": hashlib.sha256(raw_meta).hexdigest()}
ref_default = ds.groups.index("Younger") if "Younger" in ds.groups else 0
with st.sidebar:
    reference = st.selectbox("Reference group", ds.groups, index=ref_default)
    target = next(g for g in ds.groups if g != reference)
    st.caption(f"Contrast: **{target}** vs **{reference}**")

qc = qc_summary(ds)
de = differential_expression(ds, reference, target)
age = age_association(ds)
pca = sample_pca(ds)
sig = (de.fdr < fdr_threshold) & (de.log2_fold_change.abs() >= effect_threshold)

st.markdown("""
<div class="hero">
  <div class="eyebrow">COMPUTATIONAL BIOLOGY / EXPRESSION ANALYSIS</div>
  <div class="hero-title">Aging Genomics Lab</div>
  <div class="hero-sub">A reproducible research-software workbench for exploring group-level
    expression differences, age-related associations, quality control and sample structure.
    Built for transparent methods, traceable inputs and extensible bioinformatics pipelines.</div>
</div>
""", unsafe_allow_html=True)
if ds.is_simulated:
    st.markdown('<div class="warning-panel"><b>SYNTHETIC DATA / DEMONSTRATION MODE.</b> '
                'Every gene ID and expression value was generated programmatically. '
                'Charts are for software validation only—not biological discoveries.</div>',
                unsafe_allow_html=True)
else:
    st.warning("USER-IMPORTED DATA — analysis assumes correctly normalized, log2-scaled "
               "expression. Raw sequencing counts and missing values are unsupported.")

if abs(qc["age_group_correlation"]) > 0.8:
    st.caption("⚠️ This dataset has strong age–group correlation "
               f"(r={qc['age_group_correlation']:.2f}); associations with age cannot be interpreted "
               "independently from group effects without an appropriate study design/model.")

tab1, tab2, tab3, tab4 = st.tabs(["◈  Overview & QC", "↗  Differential expression",
                                   "◎  Age association", "⇩  Data & export"])

with tab1:
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Samples", f"{qc['n_samples']:,}")
    m2.metric("Genes / features", f"{qc['n_features']:,}")
    m3.metric("Significant DE features", f"{int(sig.sum()):,}")
    m4.metric("Age span", f"{qc['age_min']:.0f}–{qc['age_max']:.0f} years")
    st.write("")
    left, right = st.columns([1.45, 1])
    with left:
        st.subheader("Sample structure · principal components")
        fig = px.scatter(pca, x="pc1", y="pc2", color="group", hover_data=["sample_id", "age"],
                         color_discrete_sequence=[TEAL, PURPLE])
        fig.update_xaxes(title=f"PC1 · {pca.pc1_variance_pct.iloc[0]:.1f}% variance")
        fig.update_yaxes(title=f"PC2 · {pca.pc2_variance_pct.iloc[0]:.1f}% variance")
        st.plotly_chart(style_fig(fig), use_container_width=True)
    with right:
        st.subheader("Age distribution by cohort")
        meta_plot = ds.metadata.reset_index()
        fig = px.box(meta_plot, x="group", y="age", color="group", points="all",
                     hover_data=["sample_id"], color_discrete_sequence=[TEAL, PURPLE])
        fig.update_xaxes(title="Cohort")
        fig.update_yaxes(title="Age (years)")
        fig.update_layout(showlegend=False)
        st.plotly_chart(style_fig(fig), use_container_width=True)
    st.subheader("Input quality checks")
    c1, c2 = st.columns(2)
    with c1:
        st.dataframe(pd.DataFrame([{"Check": "Nonfinite entries", "Value": qc["nonfinite_values"]},
                                   {"Check": "Near-constant features", "Value": qc["near_constant_features"]},
                                   {"Check": "Age/group correlation", "Value": round(qc["age_group_correlation"], 3)}]),
                     hide_index=True, use_container_width=True)
    with c2:
        groups_df = ds.metadata.groupby("group", as_index=False).agg(
            samples=("age", "count"), age_mean=("age", "mean"), age_min=("age", "min"),
            age_max=("age", "max"))
        st.dataframe(groups_df, hide_index=True, use_container_width=True)

with tab2:
    st.subheader("Differential expression · volcano plot")
    st.caption(f"**{target} − {reference}** · Welch t-test · BH-FDR across all {len(de)} features. "
               "Values must already be on a normalized log2 scale; this is NOT DESeq2.")
    plot = de.copy()
    plot["neglog10_fdr"] = -np.log10(plot.fdr.clip(lower=1e-300))
    plot["classification"] = np.select(
        [sig & (plot.log2_fold_change > 0), sig & (plot.log2_fold_change < 0)],
        [f"Higher in {target}", f"Higher in {reference}"], default="Not significant")
    fig = px.scatter(plot, x="log2_fold_change", y="neglog10_fdr",
                     color="classification", hover_data=["gene_id", "fdr", "p_value"],
                     color_discrete_map={"Not significant": "#687c94",
                                         f"Higher in {target}": TEAL,
                                         f"Higher in {reference}": PURPLE}, opacity=0.82)
    if fdr_threshold > 0:
        fig.add_hline(y=-np.log10(fdr_threshold), line_dash="dash", line_color="#8c98a8")
    fig.add_vline(x=effect_threshold, line_dash="dot", line_color="#8c98a8")
    fig.add_vline(x=-effect_threshold, line_dash="dot", line_color="#8c98a8")
    fig.update_xaxes(title="Difference in mean normalized log2 expression")
    fig.update_yaxes(title="−log10(BH-adjusted p-value)")
    st.plotly_chart(style_fig(fig, 500), use_container_width=True)
    left, right = st.columns([1, 1.15])
    with left:
        st.subheader("Top 20 differential features")
        st.dataframe(de.head(20), hide_index=True, use_container_width=True,
                     column_config={"p_value": st.column_config.NumberColumn(format="%.2e"),
                                    "fdr": st.column_config.NumberColumn(format="%.2e"),
                                    "log2_fold_change": st.column_config.NumberColumn(format="%.3f")})
    with right:
        st.subheader("Expression heatmap · top 20")
        feature_ids = de.head(20).gene_id.tolist()
        mat = ds.expression.loc[feature_ids]
        z = (mat.to_numpy() - mat.to_numpy().mean(axis=1, keepdims=True)) / np.maximum(
            mat.to_numpy().std(axis=1, keepdims=True), 1e-8)
        heat = go.Figure(data=go.Heatmap(z=z, x=mat.columns.tolist(), y=feature_ids,
                                         colorscale=[[0, "#7266bd"], [0.5, "#152b47"], [1, "#65dec6"]],
                                         zmin=-2.5, zmax=2.5, colorbar=dict(title="z-score")))
        heat.update_yaxes(autorange="reversed")
        heat.update_xaxes(showticklabels=False, title="Samples (input order)")
        st.plotly_chart(style_fig(heat, 480), use_container_width=True)
    st.download_button("Download differential expression CSV", csv_bytes(de),
                       file_name="differential_expression.csv", mime="text/csv")

with tab3:
    st.subheader("Age–expression association")
    st.caption("Spearman rank correlation per feature; BH-FDR. The analysis is exploratory "
               "and does **not** adjust for batch, sex, group, cell mixture, or other covariates.")
    scope = st.selectbox("Samples to include", ["All samples"] + ds.groups)
    scope_group = None if scope == "All samples" else scope
    try:
        chosen_age = age if scope_group is None else age_association(ds, group=scope_group)
    except ValueError as exc:
        st.warning(str(exc))
        chosen_age = None
    if chosen_age is not None:
        st.metric("Age-associated features (BH-FDR < cutoff)",
                  f"{(chosen_age.fdr < fdr_threshold).sum():,}")
        left, right = st.columns([1.15, 1])
        with left:
            default_gene = chosen_age.iloc[0].gene_id
            selected_gene = st.selectbox("Inspect gene/feature", chosen_age.gene_id.tolist(),
                                         index=chosen_age.gene_id.tolist().index(default_gene))
            samples = (ds.metadata.index if scope_group is None else
                       ds.metadata.index[ds.metadata.group == scope_group])
            gene_values = pd.DataFrame({
                "sample_id": samples,
                "age": ds.metadata.loc[samples, "age"].to_numpy(),
                "group": ds.metadata.loc[samples, "group"].to_numpy(),
                "expression": ds.expression.loc[selected_gene, samples].to_numpy(),
            })
            fig = px.scatter(gene_values, x="age", y="expression", color="group",
                             hover_data=["sample_id"], color_discrete_sequence=[TEAL, PURPLE])
            fig.update_yaxes(title="Normalized log2 expression")
            fig.update_xaxes(title="Age (years)")
            st.plotly_chart(style_fig(fig), use_container_width=True)
        with right:
            st.subheader("Top age associations")
            st.dataframe(chosen_age.head(25), hide_index=True, use_container_width=True,
                         column_config={"p_value": st.column_config.NumberColumn(format="%.2e"),
                                        "fdr": st.column_config.NumberColumn(format="%.2e"),
                                        "spearman_rho": st.column_config.NumberColumn(format="%.3f")})
        st.download_button("Download age association CSV", csv_bytes(chosen_age),
                           file_name="age_association.csv", mime="text/csv")

with tab4:
    st.subheader("Data explorer & reproducibility")
    st.caption("The source and statistical assumptions are exported alongside computed results.")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Sample metadata**")
        st.dataframe(ds.metadata.reset_index(), hide_index=True, use_container_width=True)
    with c2:
        st.markdown("**Expression preview (first 15 genes, first 8 samples)**")
        st.dataframe(ds.expression.iloc[:15, :8], use_container_width=True)
    st.markdown("**Dataset manifest**")
    st.json({"source": ds.source, "data_origin": "synthetic" if ds.is_simulated else "user-imported",
             "expression_sha256": input_hashes["expression"],
             "metadata_sha256": input_hashes["metadata"],
             "n_samples": qc["n_samples"], "n_features": qc["n_features"],
             "reference": reference, "target": target,
             "version": __version__})
    st.download_button("↓  Export reproducibility package (ZIP)",
                       archive_bytes(ds, de, age, pca, qc, reference, target, input_hashes),
                       file_name="aging_genomics_analysis_results.zip", mime="application/zip",
                       type="primary")
    st.info("For real bulk RNA-seq raw counts, a count-aware workflow (DESeq2/edgeR) "
            "must be added before biological conclusions. See docs/STATISTICAL_SCOPE.md.")

st.divider()
st.caption("Aging Genomics Lab · version " + __version__ +
           " · Local demo / exploratory research software · Not for clinical use")
