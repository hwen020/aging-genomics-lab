"""Study-specific provenance and preprocessing for public human expression."""

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from .data import DataValidationError, validate_dataset
from .ingest.geo import GeoFile, acquire

ACCESSION = "GSE5086"
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5086"
PAPER = "https://doi.org/10.1371/journal.pgen.0020115"
DISCLAIMER = "https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html"
TISSUE = "Skeletal muscle taken from rectus abdominis"
OTHER_TISSUE = "Skeletal muscle taken from other anatomical region: NON-rectus abdominis"
AGE_CUTOFF = 65
MATRIX = GeoFile(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5086/matrix/GSE5086_series_matrix.txt.gz",
    "2b2d82b8a3839e8516edc1b55bca87ee84093762741ee5ed215d7838f7c1f016",
)
ANNOTATION = GeoFile(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394",
)
CARD = {
    "accession": ACCESSION, "geo_url": GEO, "paper_url": PAPER,
    "citation": "Zahn et al. (2006), PLoS Genetics 2(7):e115",
    "platform": "GPL570 Affymetrix Human Genome U133 Plus 2.0",
    "data_type": "Expression microarray",
    "source_scale": "Linear dChip PM-only normalized probe-set intensity",
    "analysis_scale": "log2 normalized probe-set intensity",
    "feature_annotation": "GPL570 GEO annotation, probe-level tests, multi-gene mappings retained",
    "access": "Public GEO download. Submitter rights may apply. No external data redistribution.",
    "access_policy_url": DISCLAIMER,
}


def parse_matrix(path: Path):
    """Read submitter sample fields and a single, rectangular GEO matrix."""
    rows = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t", strict=True)
        for key, *values in (row for row in reader if row):
            if key == "!series_matrix_table_begin":
                break
            if key in {"!Series_geo_accession", "!Sample_geo_accession", "!Sample_platform_id",
                       "!Sample_source_name_ch1", "!Sample_characteristics_ch1",
                       "!Sample_data_processing"}:
                rows[key].append(values)
        if rows["!Series_geo_accession"] != [[ACCESSION]]:
            raise DataValidationError("Expected the GSE5086 series matrix.")
        accessions = rows["!Sample_geo_accession"]
        if len(accessions) != 1:
            raise DataValidationError("Expected one GEO sample accession row.")
        samples = accessions[0]
        header = next(reader, [])
        if (not samples or any(not re.fullmatch(r"GSM\d+", x) for x in samples)
                or len(set(samples)) != len(samples) or not header or header[0] != "ID_REF"
                or len(header[1:]) != len(samples) or set(header[1:]) != set(samples)):
            raise DataValidationError("Matrix and sample metadata IDs differ or are invalid.")
        for key, fields in rows.items():
            if key != "!Series_geo_accession" and any(len(row) != len(samples) for row in fields):
                raise DataValidationError("GEO metadata row has the wrong sample count.")

        platform, tissue, processing = (
            rows[key] for key in ("!Sample_platform_id", "!Sample_source_name_ch1", "!Sample_data_processing")
        )
        if platform != [["GPL570"] * len(samples)] or len(processing) != 1 or any(
            "dchip" not in value.lower() or "pm only" not in value.lower() for value in processing[0]
        ):
            raise DataValidationError("Expected GPL570 dChip PM-only input. Counts and unknown scales are unsupported.")
        if len(tissue) != 1 or not set(tissue[0]) <= {TISSUE, OTHER_TISSUE}:
            raise DataValidationError("Unexpected or missing anatomical source.")
        fields = defaultdict(list)
        for row in rows["!Sample_characteristics_ch1"]:
            keys = {cell.partition(":")[0].strip().lower() for cell in row}
            if len(keys) == 1:
                fields[next(iter(keys))].append([cell.partition(":")[2].strip() for cell in row])
        if len(fields["age(yrs.)"]) != 1 or len(fields["sex"]) != 1:
            raise DataValidationError("Expected one age and one sex field per GEO sample.")
        age = pd.to_numeric(pd.Series(fields["age(yrs.)"][0]), errors="coerce")
        sex = fields["sex"][0]
        if not age.between(0, 120).all() or not set(sex) <= {"Female", "Male"}:
            raise DataValidationError("Invalid or missing age or sex in GEO metadata.")
        meta = pd.DataFrame({"sample_id": samples, "age": age, "sex": sex, "tissue": tissue[0],
                             "platform": "GPL570"})
        probes = []
        values = []
        for row in reader:
            if row == ["!series_matrix_table_end"]:
                break
            if len(row) != len(header) or not row[0].strip():
                raise DataValidationError("Invalid probe ID or matrix row width.")
            probes.append(row[0])
            values.append(row[1:])
        else:
            raise DataValidationError("Missing GEO matrix end marker.")
    if len(set(probes)) != len(probes):
        raise DataValidationError("Duplicate probe IDs in GEO matrix.")
    frame = pd.DataFrame(values, index=pd.Index(probes, name="gene_id"), columns=header[1:])
    frame = frame.replace({"": np.nan, "NA": np.nan, "null": np.nan})
    try:
        return frame.astype(float), meta
    except ValueError as exc:
        raise DataValidationError("Non-numeric GEO probe signal.") from exc


def prepare_study(expression, metadata):
    included = metadata.tissue.eq(TISSUE)
    cohort = metadata.loc[included].copy()
    cohort["group"] = np.where(cohort.age >= AGE_CUTOFF, "Older", "Younger")
    selected = expression.loc[:, cohort.sample_id]
    values = selected.to_numpy()
    finite_values = values[np.isfinite(values)]
    if (finite_values < 0).any() or not (finite_values % 1 != 0).any():
        raise DataValidationError("Expected nonnegative continuous dChip intensity, with fractional signals.")
    nonfinite = (~np.isfinite(values)).sum(axis=1)
    zeros = (values == 0).sum(axis=1)
    finite = nonfinite == 0
    positive = (values > 0).all(axis=1)
    logged = np.log2(selected.loc[finite & positive])
    variable = logged.std(axis=1, ddof=0) > 1e-8
    filter_log = pd.DataFrame({
        "gene_id": selected.index,
        "nonfinite_values": nonfinite,
        "zero_values": zeros,
        "exclusion_reason": np.select([~finite, ~positive], ["Nonfinite signal", "Zero intensity"], default=""),
    })
    filter_log.loc[filter_log.gene_id.isin(logged.index[~variable]), "exclusion_reason"] = "Near-constant log2 signal"
    ds = validate_dataset(logged.loc[variable].reset_index(), cohort, source=ACCESSION)
    sample_log = metadata.assign(included=included, exclusion_reason=np.where(
        included, "", "Anatomical source outside rectus abdominis"))
    return ds, sample_log, filter_log


def parse_annotation(path):
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE)
        for header in reader:
            if header and header[0] == "ID":
                break
        else:
            raise DataValidationError("Missing GPL570 annotation header.")
        columns = ["ID", "Gene symbol", "Gene ID"]
        if not set(columns) <= set(header):
            raise DataValidationError("Missing GPL570 annotation columns.")
        indices = [header.index(name) for name in columns]
        values = []
        for row in reader:
            if not row or row[0].startswith("!"):
                continue
            if len(row) <= max(indices):
                raise DataValidationError("Invalid GPL570 annotation row.")
            values.append([row[i] or None for i in indices])
    mapping = pd.DataFrame(values, columns=columns)
    mapping = mapping.rename(columns={"ID": "gene_id", "Gene symbol": "gene_symbol", "Gene ID": "entrez_id"})
    if mapping.gene_id.duplicated().any():
        raise DataValidationError("Duplicate probe IDs in GPL570 annotation.")
    return mapping


def load_study(cache=Path("data/cache"), *, offline=False):
    matrix, matrix_record = acquire(MATRIX, cache / ACCESSION, offline=offline)
    annotation, annotation_record = acquire(ANNOTATION, cache / ACCESSION, offline=offline)
    expression, metadata = parse_matrix(matrix)
    if expression.shape != (54675, 81) or metadata.tissue.eq(TISSUE).sum() != 62:
        raise DataValidationError("Unexpected GSE5086 source dimensions or rectus cohort size.")
    ds, sample_log, filter_log = prepare_study(expression, metadata)
    mapping = parse_annotation(annotation)
    return ds, sample_log, filter_log, mapping, [matrix_record, annotation_record]
