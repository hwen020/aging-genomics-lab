"""Read-only previews of versioned aggregate outputs."""

import json
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from . import __version__

SNAPSHOT = Path(__file__).with_name("portfolio_studies.json")


def render_portfolio(style_fig):
    st.title("Aging Genomics Lab")
    st.info("PUBLIC STUDY · PRECOMPUTED SNAPSHOT")
    accession = st.selectbox("GEO study", ["GSE5086", "GSE40279"])
    st.header("Expression Analysis" if accession == "GSE5086" else "Methylation Analysis")
    try:
        summary = json.loads(SNAPSHOT.read_text(encoding="utf-8"))[accession]
        if summary["analysis_version"] != __version__:
            raise ValueError("Snapshot version differs from the app version.")
    except (OSError, ValueError, KeyError):
        st.error("The versioned portfolio snapshot is unavailable.")
        return
    st.caption(f"{accession} · v{summary['analysis_version']} · Generated {summary['generated_utc']}")
    st.markdown(f"[GEO source]({summary['source']}) · [Original paper]({summary['paper']})")
    overview, analysis, results = st.tabs(["Overview", "Analysis", "Results"])
    with overview:
        a, b, c = st.columns(3)
        a.metric("Samples", f"{summary['samples']:,}")
        b.metric("Valid age tests", f"{summary['tests']:,}")
        c.metric("BH-FDR < 0.05", f"{summary['fdr_hits']:,}")
        st.code(summary["formula"], language=None)
        st.caption("Probe-wise OLS with HC3 standard errors and genome-wide Benjamini–Hochberg correction.")
        for limit in summary["limitations"]:
            st.caption(limit)
        if counts := summary.get("published_comparison"):
            st.write(f"Zahn Table S2 audit · {counts['matched_published_rows']} of {counts['published_rows']} "
                     f"published records matched, forming {counts['matched_pairs']} record–probe pairs. "
                     f"{counts['same_direction_pairs']} pairs agreed in direction.")
    with analysis:
        st.subheader("Age-association p-value distribution")
        plot = go.Figure(go.Bar(x=[(i + 0.5) / 20 for i in range(20)], y=summary["p_histogram"],
                               width=0.046, marker_color="#61D6C4"))
        plot.update_xaxes(title="Unadjusted p-value (bin midpoint)", range=[0, 1])
        plot.update_yaxes(title="Valid tests")
        st.plotly_chart(style_fig(plot), use_container_width=True)
    with results:
        st.subheader("Example age-association results")
        st.caption("Eight leading rows from the cached analysis, including non-significant rows when applicable.")
        st.dataframe(pd.DataFrame(summary["examples"]), hide_index=True, use_container_width=True)
        st.download_button("Download snapshot JSON", json.dumps(summary, indent=2),
                           f"{accession}-portfolio.json", "application/json")
