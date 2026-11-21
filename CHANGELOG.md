# Changelog

## 0.3.1 · Release preparation · 2026-11-20

- Added an aggregate-only Portfolio Demo and a repeatable snapshot builder.
- Audited both public-study models and added independent checks of estimates, robust standard errors, beta clipping and BH-FDR.
- Clarified published comparison counts for papers' rows, distinct probes and record–probe pairs.
- Replaced the EWAS tokenizer with a fixed-size beta row buffer after a Windows memory-allocation failure.
- Added citation metadata, contribution and deployment guides, a release checklist and portfolio health checking.

## 0.3.0 · P1.1 and P2 · 2026-11-20

- Added an independent Zahn Table S2 audit with published source validation, exact identifier matching, direction comparison, and report and provenance exports.
- Added a separate GSE40279 workflow with resumable beta and platform downloads, verified metadata, bounded beta blocks, probe and sample QC, and sex/plate-adjusted M-value regression.
- Added a methylation workbench with PCA, volcano, observation plots, significant CpG search and reproducibility downloads.
- Added offline coverage for invalid samples, probe filtering, model rank, chunk budgets, corrupted caches and HTTP range retries.

## 0.2.0 · P1 · 2026-11-20

- Added verified GEO acquisition and cache for the GSE5086 matrix and GPL570 annotation.
- Added rectus abdominis cohort selection, age and sex parsing, probe QC, sex-adjusted HC3 models, BH-FDR and interactive diagnostics.
- Added a Public Study workspace and shared CLI/dashboard export with CSV, Markdown, HTML charts and provenance.
- Added explicit raw-count rejection and coverage for small integer signals, data and sample validation, network failures, cache integrity and the study catalog.
- Updated the Windows launcher to reuse an installed environment and fall back to a compatible Python command during setup.

## 0.1.0

- Bundled deterministic synthetic data, normalized-expression CSV ingestion, exploratory statistics, CLI reporting and Streamlit dashboard.
