"""Streamlit workbench for the GSE40279 EWAS."""

import json
import tempfile
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd
import streamlit as st

from .methylation_analysis import FDR, run_ewas
from .methylation_reporting import export_ewas, figures
from .methylation_study import GEO, PAPER, TERMS, acquire_study

OUT = Path("output/GSE40279")


def open_report(path):
    manifest = json.loads((path / "provenance.json").read_text(encoding="utf-8"))
    if manifest["accession"] != "GSE40279":
        raise ValueError("The result directory does not contain a GSE40279 report.")
    names = ["ewas_results", "sample_pca", "top_cpg_beta", "sample_qc", "sample_metadata"]
    frames = [pd.read_csv(path / f"{name}.csv", dtype={"chromosome": str, "plate": str})
              for name in names]
    return manifest, *frames


def render_methylation(style_fig):
    st.header("Methylation Analysis")
    st.markdown(f"**GSE40279 · Whole-blood Illumina HumanMethylation450**  \n"
                f"[GEO record]({GEO}) · [Hannum et al.]({PAPER}) · [Access terms]({TERMS})")
    overview, explorer, analysis, reports = st.tabs(["Overview", "Data Explorer", "Analysis", "Results & Reports"])
    with explorer:
        st.caption("First download uses about 1.3 GB. The matrix is read in bounded CpG blocks. "
                   "Verified caches support an offline run.")
        open_existing = st.button("Open local EWAS report", disabled=not (OUT / "provenance.json").is_file())
        rerun = st.button("Verify cache and rerun EWAS")
        download = st.button("Download GEO data and run EWAS")
        if open_existing or rerun or download:
            try:
                with st.spinner("Opening report or verifying files and fitting EWAS models…"):
                    if rerun or download:
                        paths, records = acquire_study(offline=rerun)
                        export_ewas(run_ewas(paths), records, OUT)
                    st.session_state["gse40279"] = open_report(OUT)
                    st.session_state.pop("gse40279_zip", None)
            except (OSError, ValueError, EOFError, KeyError) as exc:
                st.error(str(exc))
    state = st.session_state.get("gse40279")
    if state is None:
        with overview:
            left, right = st.columns(2)
            left.metric("Submitted samples", 656)
            right.metric("Submitted scale", "Beta fractions")
            st.info("Open Data Explorer to load or generate a local report.")
        return

    manifest, results, pca, betas, sample_qc, meta = state
    qc = manifest["qc"]
    with overview:
        a, b, c = st.columns(3)
        a.metric("Analyzed samples", f"{qc['samples']:,}")
        b.metric("Tested CpGs", f"{qc['tested_probes']:,}")
        c.metric("FDR < 0.05", f"{qc['fdr_hits']:,}")
        st.write("How does whole-blood methylation vary with age after accounting for recorded sex and plate?")
        st.caption(f"Age range {qc['age_range'][0]}–{qc['age_range'][1]} years. Input and plots use beta fractions. "
                   "Inference uses log₂ beta odds per decade.")
    with explorer:
        st.subheader("Samples and phenotypes")
        st.dataframe(meta, hide_index=True)
        st.subheader("Probe exclusions")
        st.table(pd.Series(qc["excluded_probes"], name="Probe count"))
        st.subheader("Sample QC")
        st.dataframe(sample_qc, hide_index=True)
        with st.expander("Verified inputs and SHA-256"):
            st.json(manifest["downloaded_artifacts"])
    with analysis:
        st.subheader("Sample PCA")
        st.caption(f"PCA uses {manifest['parameters']['pca_probes']:,} deterministically selected eligible beta probes.")
        probe = st.selectbox("CpG to plot", betas.cpg_id.unique())
        charts = figures(results, pca, betas, probe=probe)
        st.plotly_chart(style_fig(charts["pca"]), use_container_width=True)
        st.subheader("Continuous age association")
        st.caption("Volcano displays the strongest 1,000 tests plus a deterministic subset. All tests are exported.")
        st.plotly_chart(style_fig(charts["volcano"]), use_container_width=True)
        st.caption("The scatter shows observed beta fractions. Model effects use M-values and adjust for sex and plate.")
        st.plotly_chart(style_fig(charts["age_beta"]), use_container_width=True)
    with reports:
        st.subheader("EWAS results")
        query = st.text_input("Search CpG or gene symbol").strip()
        only_hits = st.checkbox("Significant CpGs only", value=True)
        subset = results.loc[results.fdr < FDR] if only_hits else results
        if query:
            subset = subset.loc[subset.cpg_id.str.contains(query, case=False, regex=False)
                                | subset.gene_symbols.fillna("").str.contains(query, case=False, regex=False)]
        st.caption(f"{len(subset):,} matching rows. Showing up to 500.")
        st.dataframe(subset.head(500), hide_index=True)
        st.download_button("Download filtered CSV", subset.to_csv(index=False), "GSE40279-filtered.csv", "text/csv")
        st.markdown((OUT / "report.md").read_text(encoding="utf-8"))
        st.download_button("Download Markdown report", (OUT / "report.md").read_bytes(), "GSE40279-report.md")
        if st.button("Build EWAS reproducibility ZIP"):
            with tempfile.TemporaryFile() as buffer:
                with ZipFile(buffer, "w", ZIP_DEFLATED) as archive:
                    for item in sorted(OUT.iterdir()):
                        if item.suffix in {".csv", ".html", ".md", ".json"} and item.is_file():
                            archive.write(item, item.name)
                buffer.seek(0)
                st.session_state["gse40279_zip"] = buffer.read()
        if package := st.session_state.get("gse40279_zip"):
            st.download_button("Download EWAS ZIP", package, "GSE40279-report.zip", "application/zip")
