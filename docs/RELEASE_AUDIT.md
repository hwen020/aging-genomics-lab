# v0.3.1 scientific audit

## Alignment and quality control

GSE5086 joins the deposited probe matrix to its GSM metadata by unique sample identifiers and orders the expression columns against the selected cohort. This leaves 62 rectus abdominis samples and 54,615 tested GPL570 probe sets. The 60 exclusions for zero intensity occur before log2. The age slope is expression change per decade conditional on sex. The group contrast uses Older minus Younger conditional on sex. Age defines that contrast, and GEO provides no batch labels.

GSE40279 links GSM to submitter IDs and the official sample key, validates matrix column uniqueness and aligns each beta block to metadata order. QC is complete-case by CpG and checks beta range, annotation, nearby annotated SNPs and minimum variability. Plate has nine observed levels. The model uses an intercept, centered age per decade, sex and eight plate indicators with plate 1 as reference. The verified matrix has rank 11 and residual df 645. The implementation clips only for log2-odds inference at 1e-6 and 1−1e-6. Plots use submitted beta values.

## Independent inference checks

The release tests independently calculate OLS via normal equations, construct the scalar HC3 sandwich from individual observations and evaluate the t-tail with the incomplete beta function. Fixtures check signed age and group effects, beta endpoints, plate coding, rank failure and an unsorted BH example with ties. The release audit independently matched six expression estimates, six methylation estimates, their errors and p-values, and BH over all valid tests in both expression families and the EWAS family.

The prior beta parser twice reported a Windows allocation error on the full cached matrix. The release reader validates tab-delimited records and uses a fixed-size numeric buffer. With this fix the entire EWAS completed at the default block size. Compared with v0.3.0, both studies had numerically equivalent estimates, PCA and QC and identical sample/probe sets and FDR decisions. Sign alignment was allowed for PCA axes, and numeric comparisons used tolerance for float parsing and BLAS differences.

## Published comparison

Zahn's 250 findings are represented by published rows. The comparison links 113 distinct published records to 219 distinct GPL570 probes and yields 219 record–probe pairs, with 196 same-direction pairs. Exact probe IDs take precedence and single exact gene symbols may yield multiple probe links. This is a descriptive comparison using overlapping source samples. The analysis models, thresholds and matching rules retain their original definitions.

## Remaining limitations

Whole-blood cell composition, other unmeasured confounders and incomplete detection, bead-count and cross-reactivity QC may affect EWAS associations. Ethnicity is nested within plate, and recruitment source overlaps with plate. GPL13534 gene annotations describe genomic context and can be incomplete or ambiguous. The cross-sectional studies support exploratory associations. The probe-level GSE5086 models do not replicate the paper's gene aggregation and permutation decision rule.
