# Statistical scope and limitations

## P1 public microarray analysis

GSE5086 uses 62 rectus abdominis samples selected by the submitted anatomical source. The 19 excluded samples remain in the sample log. GSM identifiers match between the metadata and expression table. The source is dChip PM-only normalized linear probe intensity. Positive finite values are transformed once with log2. A probe is excluded if any retained sample has a nonfinite or zero intensity or if its log2 population standard deviation is at most 1e-8. No values are imputed. Negative intensities, integer-only signals, unknown processing, incomplete age/sex metadata and sample mismatches terminate the run. No extra normalization or batch correction is applied.

Each model tests one coefficient per GPL570 probe set. The age model uses centered age per decade and the group model uses age ≥65 versus age <65. Both include an intercept and submitter-recorded sex with female as the reference. Ordinary least squares estimates use HC3 standard errors, two-sided t tests and 95% intervals with residual degrees of freedom. A rank-deficient design or leverage approaching one stops the analysis. A probe with nonfinite HC3 error or error ≤1e-12 has no p-value or FDR. Benjamini–Hochberg operates separately on valid p-values from each model. Report defaults are FDR <0.05 and, for the group model, absolute log2 effect ≥0.5. Heatmap selection uses the top 20 valid group results by FDR, with within-probe z scores for display.

Centered sample PCA and observed expression-by-age scatter are descriptive. QC flags sample median log2 signals with absolute median/MAD robust z score >3.5 using a 1.4826 scaling factor. The score is unavailable at zero MAD. Flags are retained in the cohort. Batch labels were unavailable. Age and age group describe the same exposure, and unmeasured covariates may confound estimates. A probe may map to several or no genes. The cutoff, cohort and models differ from the original publication and support exploratory associations only. [Source and access details](PUBLIC_STUDIES.md)

## Synthetic demo and normalized CSV workflow

**Group comparison**: Welch's two-sample t-test for each gene on already normalized, log2-scale continuous expression, followed by Benjamini–Hochberg false discovery rate adjustment. The `log2_fold_change` column is the arithmetic difference in group means on the input log2 scale, not a count-model estimate.

**Age association**: per-feature Spearman rank correlation of age with input expression, followed by BH FDR. The p-value is exploratory and is *not* adjusted for group, batch, sex, ancestry, inflammation, comorbidities, tissue composition, repeated measures or other confounding.

**PCA**: gene-centered sample-level singular-value decomposition; scores and explained variance percentages. It is unsupervised, but any apparent separation may reflect technical artifacts.

**QC**: finite values, duplicate/mismatched sample ID validation, near-constant feature count, group age statistics and age–group correlation.

## Not implemented / don't claim

- DESeq2/edgeR negative binomial analysis of raw RNA-seq counts.
- Random effects for donors, longitudinal repeated measures and covariate adjustment in the legacy CSV workflow.
- Batch-effect correction, surrogate variable analysis, quantile normalization, RMA or RNA-seq read alignment.
- DNA methylation EWAS or epigenetic clocks, genetic variant calling, GWAS, TWAS.
- Gene enrichment, regulatory motif learning, validated disease biomarker discovery.
- Any clinical diagnosis, personalized disease risk prediction or causal discovery.

## Highest-risk caveat

The synthetic demo intentionally has near-separate age ranges across the Younger and Older cohorts. Age and group are therefore **strongly confounded**. Correlation with age must not be described as an independent aging mechanism, including when Spearman p-values appear small. This project is primarily evidence of data-engineering, statistics implementation, testing and reproducibility.

## Future count workflow

Use a real public count dataset (e.g. GSE164471), document sample selection and donor-level covariates, inspect distributions and independent biological replicates, fit a DESeq2 design accounting for batch/sex when appropriate, report design matrix and contrasts, control FDR, and independently reproduce published patterns with cautious language. Include objective negative-control checks and uncertainty reporting.
