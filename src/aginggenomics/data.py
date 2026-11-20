"""Strict ingestion and validation for pre-normalized log2 gene expression."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import IO

import numpy as np
import pandas as pd


class DataValidationError(ValueError):
    """Raised for unsafe/inconsistent experiment input."""


@dataclass(frozen=True)
class Dataset:
    expression: pd.DataFrame  # gene_id index, sample_id columns; values log2 normalized
    metadata: pd.DataFrame    # sample_id index; group, age columns
    source: str = "user-imported"
    is_simulated: bool = False

    @property
    def groups(self) -> list[str]:
        return sorted(self.metadata["group"].astype(str).unique().tolist())


def load_dataset(expression_file: str | Path | IO, metadata_file: str | Path | IO,
                 *, source: str = "user-imported", simulated: bool = False) -> Dataset:
    """Read two CSVs with matching sample IDs; never silently impute or normalize."""
    try:
        expr = pd.read_csv(expression_file, dtype={"gene_id": str})
        meta = pd.read_csv(metadata_file, dtype={"sample_id": str, "group": str})
    except Exception as exc:
        raise DataValidationError(f"Unable to parse CSV inputs: {exc}") from exc
    return validate_dataset(expr, meta, source=source, simulated=simulated)


def validate_dataset(expression: pd.DataFrame, metadata: pd.DataFrame,
                     *, source: str = "user-imported", simulated: bool = False) -> Dataset:
    expression = expression.copy()
    metadata = metadata.copy()
    if "gene_id" not in expression.columns:
        raise DataValidationError("Expression CSV must have a 'gene_id' column.")
    needed = {"sample_id", "group", "age"}
    if not needed.issubset(metadata.columns):
        raise DataValidationError(f"Metadata CSV must have columns: {', '.join(sorted(needed))}.")
    if expression.shape[1] < 5:
        raise DataValidationError("Expression CSV must contain at least four sample columns.")
    if len(expression) < 2:
        raise DataValidationError("At least two genes/features are required.")

    expression["gene_id"] = expression["gene_id"].astype(str).str.strip()
    metadata["sample_id"] = metadata["sample_id"].astype(str).str.strip()
    metadata["group"] = metadata["group"].astype(str).str.strip()
    if expression["gene_id"].eq("").any() or expression["gene_id"].duplicated().any():
        raise DataValidationError("Gene IDs must be nonempty and unique.")
    if metadata["sample_id"].eq("").any() or metadata["sample_id"].duplicated().any():
        raise DataValidationError("Sample IDs must be nonempty and unique in metadata.")
    if metadata["group"].eq("").any() or (metadata["group"].str.lower() == "nan").any():
        raise DataValidationError("Group values must be nonempty.")

    samples = [str(c).strip() for c in expression.columns if c != "gene_id"]
    if len(set(samples)) != len(samples) or any(not sample for sample in samples):
        raise DataValidationError("Expression sample columns must be nonempty and unique.")
    metadata_samples = metadata["sample_id"].tolist()
    if set(samples) != set(metadata_samples):
        missing_meta = sorted(set(samples) - set(metadata_samples))
        missing_expr = sorted(set(metadata_samples) - set(samples))
        raise DataValidationError(
            f"Sample IDs differ between files. Missing metadata: {missing_meta[:8]}; "
            f"missing expression: {missing_expr[:8]}"
        )

    try:
        data = expression.drop(columns="gene_id").apply(pd.to_numeric, errors="raise")
        values = data.to_numpy(dtype=float)
    except (ValueError, TypeError) as exc:
        raise DataValidationError("Expression contains nonnumeric values.") from exc
    if not np.isfinite(values).all():
        raise DataValidationError("Expression contains missing, NaN, or infinite values. No imputation is performed.")
    if np.nanpercentile(values, 99) > 40:
        raise DataValidationError(
            "Values look too large for normalized log2 expression; raw counts are NOT supported. "
            "Use an appropriate RNA-seq count model (e.g. DESeq2) instead."
        )
    if not simulated and (values >= 0).all() and np.equal(values, np.floor(values)).all():
        raise DataValidationError(
            "Entirely nonnegative integer expression is ambiguous with raw counts. "
            "Provide continuous normalized log2 expression or use a count-aware model."
        )
    data.index = expression["gene_id"].values
    data.index.name = "gene_id"
    data.columns = samples

    try:
        metadata["age"] = pd.to_numeric(metadata["age"], errors="raise")
    except (TypeError, ValueError) as exc:
        raise DataValidationError("Age must be numeric.") from exc
    if not np.isfinite(metadata["age"].to_numpy(dtype=float)).all():
        raise DataValidationError("Age values must be finite.")
    if not metadata["age"].between(0, 120).all():
        raise DataValidationError("Age must be in 0–120 years for this demonstration.")
    if metadata["age"].nunique() < 3:
        raise DataValidationError("At least three distinct ages are needed for age-association analysis.")

    metadata = metadata.set_index("sample_id", drop=True).loc[samples]
    metadata.index.name = "sample_id"
    counts = metadata["group"].value_counts()
    if len(counts) != 2:
        raise DataValidationError("Exactly two groups are required for the differential-expression comparison.")
    if (counts < 3).any():
        raise DataValidationError("Each group requires at least 3 biological samples.")
    return Dataset(data, metadata, source=source, is_simulated=simulated)
