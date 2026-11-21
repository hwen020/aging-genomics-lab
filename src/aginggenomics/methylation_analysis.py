"""Blockwise continuous-age EWAS using logit-transformed methylation beta values."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.stats import t

from .analysis import benjamini_hochberg
from .data import DataValidationError
from .methylation_study import AUTOSOMES, PROBE_COUNT, beta_chunks, parse_phenotypes, parse_platform

BETA_EPSILON = 1e-6
MIN_BETA_SD = 0.005
FDR = 0.05
PCA_MODULUS = 64
TOP_PROBES = 20
MODEL_FORMULA = "M ~ centered_age_per_decade + sex + C(plate)"
CAVEATS = [
    "Whole-blood cell proportions are unavailable. Age-related changes in cell mixture may confound associations.",
    "The model adjusts for submitted sex and plate. Ethnicity is nested within plate in this cohort, and recruitment source and plate overlap strongly.",
    "Processed beta values support range, missingness, variability and annotation QC. Detection p-values, bead counts and cross-reactivity masks were unavailable to this workflow.",
    "Results reflect a cross-sectional association study. They do not estimate biological age or establish causality.",
]


@dataclass
class EwasRun:
    metadata: pd.DataFrame
    results: pd.DataFrame
    probe_qc: pd.DataFrame
    sample_qc: pd.DataFrame
    pca: pd.DataFrame
    example_betas: dict
    parameters: dict


def age_design(metadata):
    age = metadata.age.to_numpy(dtype=float)
    plate = pd.get_dummies(metadata.plate, drop_first=True, dtype=float)
    design = np.column_stack((np.ones(len(age)), (age - age.mean()) / 10, metadata.sex.eq("M"), plate))
    df = len(age) - design.shape[1]
    if df < 3 or np.linalg.matrix_rank(design) < design.shape[1]:
        raise DataValidationError("Rank-deficient age, sex and plate design or insufficient residual degrees of freedom.")
    inverse = np.linalg.pinv(design)
    leverage = np.einsum("ij,ji->i", design, inverse)
    if (leverage >= 1 - 1e-10).any():
        raise DataValidationError("Sex/plate design contains sample leverage too high for HC3.")
    return design, inverse, inverse[1] / (1 - leverage), df


def fit_block(beta, model):
    """HC3 t inference on log2(beta / (1-beta)); beta itself stays unmodified."""
    design, inverse, weights, df = model
    clipped = np.clip(beta, BETA_EPSILON, 1 - BETA_EPSILON)
    outcome = np.log2(clipped / (1 - clipped)).T
    coefficients = inverse @ outcome
    residual = outcome - design @ coefficients
    se = np.linalg.norm(weights[:, None] * residual, axis=0)
    valid = np.isfinite(se) & (se > 1e-12)
    p_values = np.full(len(se), np.nan)
    p_values[valid] = 2 * t.sf(np.abs(coefficients[1, valid] / se[valid]), df)
    se[~valid] = np.nan
    return coefficients[1], se, p_values


def qc_block(beta, annotations, ids):
    raw_missing = np.isnan(beta)
    outside = ~raw_missing & (~np.isfinite(beta) | (beta < 0) | (beta > 1))
    missing = raw_missing.sum(axis=1)
    invalid = outside.sum(axis=1)
    annotated = (annotations.build.eq("37") & annotations.position.notna()).to_numpy()
    autosomal = annotations.chromosome.isin(AUTOSOMES).to_numpy()
    nearby_snp = annotations.nearby_snp.eq(True).to_numpy()
    with np.errstate(invalid="ignore"):
        spread = np.std(beta, axis=1)
    reason = np.select(
        [~ids.str.fullmatch(r"cg\d{8}").to_numpy(dtype=bool), invalid > 0, missing > 0,
         ~annotated, ~autosomal, nearby_snp, spread < MIN_BETA_SD],
        ["Non-CpG identifier", "Out-of-range or infinite beta", "Missing beta",
         "Missing hg19 annotation", "Non-autosomal probe", "Nearby annotated SNP",
         "Low beta variability"],
        default="",
    )
    return raw_missing, outside, missing, invalid, reason


def run_ewas(paths, *, chunk_rows=512):
    matrix, soft, key, platform = paths
    metadata = parse_phenotypes(soft, key)
    return analyze_matrix(matrix, metadata, parse_platform(platform), chunk_rows=chunk_rows,
                          expected_rows=PROBE_COUNT)


def analyze_matrix(matrix, metadata, annotations, *, chunk_rows=512, expected_rows):
    model = age_design(metadata)
    df = model[-1]
    n = len(metadata)
    seen = set()
    gram = np.zeros((n, n))
    sample_totals = np.zeros((3, n))
    tests, qc_frames, pca_blocks = [], [], []
    top = []
    pca_count = 0
    pca_pending = 0
    for ids, beta in beta_chunks(matrix, metadata, chunk_rows=chunk_rows):
        if ids.isna().any() or ids.duplicated().any() or seen.intersection(ids):
            raise DataValidationError("Missing or duplicate probe identifier in the beta matrix.")
        seen.update(ids)
        if len(seen) > expected_rows:
            raise DataValidationError("Beta matrix has more probes than expected.")
        anno = annotations.reindex(ids)
        missing_cells, invalid_cells, nm, ni, reason = qc_block(beta, anno, ids)
        sample_totals += np.vstack((missing_cells.sum(axis=0), invalid_cells.sum(axis=0),
                                    np.where(missing_cells | invalid_cells, 0, beta).sum(axis=0)))
        qc_frames.append(pd.DataFrame({"cpg_id": ids.to_numpy(), "missing_values": nm,
                                       "invalid_values": ni, "exclusion_reason": reason}))
        keep = reason == ""
        if not keep.any():
            continue
        kept_ids = ids.to_numpy()[keep]
        kept = beta[keep]
        effect, se, p = fit_block(kept, model)
        critical = t.ppf(0.975, df)
        tests.append(pd.DataFrame({"cpg_id": kept_ids, "m_change_per_decade": effect,
                                   "standard_error": se, "ci95_lower": effect - critical * se,
                                   "ci95_upper": effect + critical * se, "p_value": p}))
        eligible = [i for i in np.argsort(p, kind="stable")[:TOP_PROBES] if np.isfinite(p[i])]
        top.extend((float(p[i]), kept_ids[i], kept[i].copy()) for i in eligible)
        top = sorted(top, key=lambda item: item[:2])[:TOP_PROBES]
        pca_mask = np.fromiter((int(identifier[2:]) % PCA_MODULUS == 0 for identifier in kept_ids),
                               dtype=bool, count=len(kept_ids))
        subset = kept[pca_mask]
        if len(subset):
            pca_blocks.append(subset - subset.mean(axis=1, keepdims=True))
            pca_count += len(subset)
            pca_pending += len(subset)
        if pca_pending >= 256:
            block = np.vstack(pca_blocks)
            gram += block.T @ block
            pca_blocks.clear()
            pca_pending = 0

    if len(seen) != expected_rows or not tests or pca_count < 2:
        raise DataValidationError("Unexpected probe count, or insufficient eligible probes for EWAS/PCA.")
    if pca_blocks:
        block = np.vstack(pca_blocks)
        gram += block.T @ block
    eigenvalues, vectors = eigh(gram, subset_by_index=(n - 2, n - 1))
    scores = vectors[:, ::-1] * np.sqrt(np.maximum(eigenvalues[::-1], 0))
    for column in scores.T:
        column *= np.sign(column[np.argmax(np.abs(column))]) or 1
    pc1_pct, pc2_pct = 100 * np.maximum(eigenvalues[::-1], 0) / np.trace(gram)
    pca = metadata[["sample_id", "age", "sex", "plate"]].assign(
        pc1=scores[:, 0], pc2=scores[:, 1], pc1_variance_pct=pc1_pct, pc2_variance_pct=pc2_pct)
    results = pd.concat(tests, ignore_index=True)
    valid = results.p_value.notna()
    results.loc[valid, "fdr"] = benjamini_hochberg(results.loc[valid, "p_value"].to_numpy())
    results = results.join(annotations[["chromosome", "position", "gene_symbols", "island_relation"]],
                           on="cpg_id").sort_values(["fdr", "p_value", "cpg_id"], ignore_index=True)
    missing, invalid, totals = sample_totals
    observed = expected_rows - missing - invalid
    means = np.divide(totals, observed, out=np.full(n, np.nan), where=observed > 0)
    median = np.nanmedian(means)
    mad = np.nanmedian(np.abs(means - median))
    robust_z = (means - median) / (1.4826 * mad) if mad > 0 else np.full(n, np.nan)
    sample_qc = metadata[["sample_id"]].assign(missing_values=missing.astype(int),
        out_of_range_values=invalid.astype(int), mean_beta=means, mean_beta_robust_z=robust_z,
        flag_outlier=np.abs(robust_z) > 3.5, flag_missing=(missing + invalid) / expected_rows > 0.05)
    parameters = {"model_formula": MODEL_FORMULA, "residual_df": df, "fdr": FDR,
                  "beta_clip_epsilon": BETA_EPSILON, "min_beta_sd": MIN_BETA_SD,
                  "chunk_rows": chunk_rows, "pca_probes": pca_count, "pca_probe_id_modulus": PCA_MODULUS,
                  "sample_mean_beta_flag_abs_robust_z": 3.5, "sample_missingness_flag_fraction": 0.05,
                  "design_condition_number": float(np.linalg.cond(model[0])), "random_seed": None}
    return EwasRun(metadata, results, pd.concat(qc_frames, ignore_index=True), sample_qc, pca,
                   {identifier: beta for _, identifier, beta in top}, parameters)
