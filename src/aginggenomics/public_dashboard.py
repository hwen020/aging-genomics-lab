"""Streamlit view of the verified public microarray workflow."""

from io import BytesIO
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd
import streamlit as st

from .data import DataValidationError
from .public_analysis import LIMITATIONS, MODEL, analyze_study, build_export, make_figures
from .public_studies import ACCESSION, ANNOTATION, CARD, MATRIX, load_study


def render_public_study(style_fig):
    st.title("Aging Genomics Lab")
    st.success("PUBLIC STUDY · Human skeletal muscle expression")
    selection = st.selectbox("GEO study", [ACCESSION, "GSE164471"])
    if selection == "GSE164471":
        st.info("GSE164471 · Skeletal Muscle Transcriptome in Healthy Aging. Gene-level RNA-seq counts "
                "require an independent count-aware workflow such as DESeq2.")
        st.markdown("[GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164471) · "
                    "[Original paper](https://doi.org/10.1038/s41467-021-22168-2)")
        return
    st.markdown(f"**{CARD['platform']}** · {CARD['citation']}  \n"
                f"[GEO record]({CARD['geo_url']}) · [Paper]({CARD['paper_url']}) · "
                f"[Access terms]({CARD['access_policy_url']})")
    overview, explorer, analysis, reports = st.tabs(["Overview", "Data Explorer", "Analysis", "Results & Reports"])
    with explorer:
        cache = Path("data/cache")
        local = all((cache / ACCESSION / artifact.filename).is_file() for artifact in (MATRIX, ANNOTATION))
        status = st.empty()
        load_cache = st.button("Verify and analyze cached files", disabled=not local)
        download = st.button("Download from GEO and analyze")
        if load_cache or download:
            try:
                with st.spinner("Verifying files and fitting probe-level models…"):
                    study = load_study(cache, offline=load_cache)
                    tables = analyze_study(study)
                st.session_state["gse5086_run"] = (study, tables)
                st.session_state.pop("gse5086_export", None)
                local = True
            except (OSError, DataValidationError) as exc:
                st.error(str(exc))
                if "gse5086_run" in st.session_state:
                    st.caption("The previous verified analysis remains displayed.")
        status.caption("Local files present. SHA-256 is verified on load." if local else "Official GEO download required.")
    if "gse5086_run" not in st.session_state:
        with overview:
            a, b, c, d = st.columns(4)
            a.metric("Submitted samples", 81)
            b.metric("Submitted probe sets", "54,675")
            c.metric("Data type", "Microarray")
            d.metric("Analysis", "Ready to load")
            st.info("Open Data Explorer to load the GEO data. The analysis is restricted to rectus abdominis.")
        return

    study, tables = st.session_state["gse5086_run"]
    ds, samples, filters, _, files = study
    group, age = tables["differential_expression"], tables["age_association"]
    with overview:
        a, b, c, d = st.columns(4)
        a.metric("Analyzed samples", len(ds.metadata))
        b.metric("Retained probe sets", f"{len(ds.expression):,}")
        c.metric("Data type", "Microarray")
        d.metric("Analysis", "Complete")
        st.markdown("**Research question** · How does normalized expression in rectus abdominis vary with age after accounting for sex?")
        st.caption(CARD["source_scale"] + " · Converted once to log2 after filtering.")
        st.caption("Age groups use a 65-year cutoff. Group differences and continuous-age associations share the same age exposure.")
        st.table(ds.metadata.groupby("group").age.agg(["count", "min", "median", "max"]))
    with explorer:
        st.success("Source files passed SHA-256 verification.")
        with st.expander("Source hashes and download times"):
            st.json(files)
        st.subheader("Cohort and sample exclusions")
        st.dataframe(samples, hide_index=True, use_container_width=True)
        st.subheader("Expression preview · log2 signal")
        st.dataframe(ds.expression.iloc[:15, :10], use_container_width=True)
        st.subheader("Probe and sample QC")
        st.caption(f"Before probe filtering in the selected cohort: {filters.nonfinite_values.sum():,} nonfinite "
                   f"and {filters.zero_values.sum():,} zero intensity values.")
        st.dataframe(filters.exclusion_reason.replace("", "Retained").value_counts().rename("Probe sets"))
        st.caption("Sample median flags use |robust z| > 3.5. Flagged samples remain in the model.")
        st.dataframe(tables["sample_qc"], hide_index=True, use_container_width=True)
    with reports:
        fdr = st.slider("FDR cutoff", 0.01, 0.20, 0.05, 0.01)
        effect = st.slider("Minimum absolute group log2 effect", 0.0, 2.0, 0.5, 0.1)
        view = st.radio("Result table", ["Age association", "Differential expression"], horizontal=True)
        chosen = age if view == "Age association" else group
        search = st.text_input("Find probe ID or gene symbol").strip()
        hits = chosen.fdr < fdr
        if view == "Differential expression":
            hits &= chosen.log2_fold_change.abs() >= effect
        st.metric("Probes passing thresholds", f"{int(hits.sum()):,}")
        subset = chosen[hits] if st.checkbox("Show only probes passing thresholds") else chosen
        if search:
            subset = subset[subset.gene_id.str.contains(search, case=False, regex=False)
                            | subset.gene_symbol.fillna("").str.contains(search, case=False, regex=False)]
        st.caption(f"{len(subset):,} matching rows · showing up to 500")
        st.dataframe(subset.head(500), hide_index=True, use_container_width=True)
        st.download_button("Download filtered CSV", subset.to_csv(index=False), "filtered_results.csv", "text/csv")
        st.write(MODEL + ". Benjamini–Hochberg FDR is computed separately for each model.")
        for limitation in LIMITATIONS:
            st.caption(limitation)
    with analysis:
        st.subheader("Sample structure and age distribution")
        first, second = st.columns(2)
        age_probes = age.dropna(subset=["p_value"]).gene_id
        probe_options = age_probes.head(25).tolist()
        requested_probe = st.text_input("Plot a specific probe ID").strip()
        if requested_probe in ds.expression.index:
            probe_options = list(dict.fromkeys([requested_probe, *probe_options]))
        elif requested_probe:
            st.warning("That probe is absent from the retained expression matrix.")
        probe = st.selectbox("Age-expression probe", probe_options) if probe_options else None
        figures = make_figures(ds, tables, fdr=fdr, effect=effect, probe=probe)
        first.plotly_chart(style_fig(figures["pca"]), use_container_width=True)
        second.plotly_chart(style_fig(figures["age_distribution"]), use_container_width=True)
        st.subheader("Covariate balance")
        st.dataframe(pd.crosstab(ds.metadata.group, ds.metadata.sex), use_container_width=True)
        st.caption("Sex is included in both models. The series matrix provides no batch labels.")
        st.subheader("Differential expression")
        st.plotly_chart(style_fig(figures["volcano"]), use_container_width=True)
        st.caption("Heatmap shows the top 20 group-model probes by FDR, standardized within each row.")
        st.plotly_chart(style_fig(figures["heatmap"]), use_container_width=True)
        if "age_expression" in figures:
            st.subheader("Expression by age")
            st.caption("The scatter shows observed values. The results table gives sex-adjusted effects per decade.")
            st.plotly_chart(style_fig(figures["age_expression"]), use_container_width=True)
    with reports:
        options = (fdr, effect, probe)
        if st.button("Build reproducibility package"):
            with st.spinner("Building reports, data tables and interactive chart downloads…"):
                artifacts, _ = build_export(study, tables, fdr=fdr, effect=effect, figures=figures, age_probe=probe)
                buffer = BytesIO()
                with ZipFile(buffer, "w", ZIP_DEFLATED) as archive:
                    for filename, contents in artifacts.items():
                        archive.writestr(filename, contents)
                st.session_state["gse5086_export"] = (options, buffer.getvalue())
        saved = st.session_state.get("gse5086_export")
        if saved and saved[0] == options:
            st.download_button("Download reproducibility ZIP", saved[1], "GSE5086-report.zip", "application/zip")
