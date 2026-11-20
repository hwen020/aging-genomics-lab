# Coding agent guardrails

You are extending an **existing, runnable v0.1.0** genomics research-software repository. Do **not** recreate the repository, replace the current app, delete the synthetic demo or break the CLI. Preserve current data schemas, results and tests unless a migration is written and validated. Prefer small, testable patches.

## Safety and research integrity

1. Never call synthetic demo samples real patients or real biological findings.
2. Keep raw RNA-seq counts separate from normalized log2 expression; use DESeq2/edgeR for raw counts, with proper replicates and design.
3. For real public data, preserve accession, license, processing provenance, source identifiers, checksums and citation.
4. Never insert private human genetic or clinical data, API keys, identifiable subjects or secrets into git.
5. Flag age/group confounding and unadjusted tests. Never make causal or clinical claims.
6. No irreversible changes to original data; ingest into a cache directory excluded from git.
7. Run `pytest`, CLI smoke test, and check Streamlit start before declaring done.
8. Document every new analysis method, parameters and failure modes; add regression tests.

## Priority

Read `CODEX_NEXT_STEPS.md`. Implement only the current phase unless the user changes scope. At completion report changed files, commands run, test totals, screenshots if practical, limitations and next milestone.
