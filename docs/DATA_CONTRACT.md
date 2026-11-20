# Data contract v0.1

## Expression matrix

- UTF-8 comma-separated `expression.csv`.
- First feature identifier column must be `gene_id` (this can be a transcript/probe ID but should be documented).
- Remaining header columns: unique sample IDs.
- Numeric **normalized log2** measurements, one row per feature; no missing or infinite values.
- Raw RNA-seq counts, UMI counts, untransformed TPM/FPKM and methylation beta values are **not supported** by this exploratory DE model.
- Sample ordering is taken from the expression matrix and metadata is explicitly reindexed to match it.

## Metadata

- UTF-8 CSV with at least `sample_id`, `group`, `age`.
- Exactly two groups, each with >=3 biological samples.
- Age in years, 0–120, numeric; at least three distinct ages for the default age-association UI.
- Extra fields such as `sex`, `batch`, `donor_id`, `tissue`, `condition`, `platform` are retained but are **not modeled** at this stage.
- Duplicate features or sample IDs, mismatches, nonnumeric expression or missing data cause explicit errors.

## Public GEO ingestion in P1

When converting a real GEO study, record: GEO accession, paper DOI, download time, file SHA-256, original sample accession, instrument/platform, processing method, reference genome/annotation, sex, tissue, phenotype, missing values, batch and exclusion reasons. **Never interpret an experimental timepoint or technical replicate as an independent human participant.**

GSE5086 is routed through a pinned study-specific parser. The source uses GPL570 probe IDs. Exports keep the existing `gene_id` result key and add platform `gene_symbol` and `entrez_id` annotations. Tables preserve GSM identifiers and the inclusion and exclusion log. The source artifacts stay under `data/cache/`.

A public RNA-seq raw-count matrix must be routed to the planned count-aware DESeq2/edgeR backend instead of being silently log2-transformed for Welch tests. The current cutoff check for huge matrix values is merely a heuristic, **not** a sufficient raw-count detector. The UI requires an explicit normalized-log2-scale confirmation for imported datasets; the CLI requires `--ack-normalized-log2` when the data are not marked synthetic.
