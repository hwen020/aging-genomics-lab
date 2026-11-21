"""GSE40279 source artifacts, phenotype matching and bounded matrix ingestion."""

import csv
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .data import DataValidationError
from .ingest.geo import GeoFile, acquire
from .ingest.resumable import acquire_large

ACCESSION = "GSE40279"
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279"
PAPER = "https://doi.org/10.1016/j.molcel.2012.10.016"
TERMS = "https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html"
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/suppl/"
MATRIX_BYTES = 1_223_258_643
MATRIX = GeoFile(BASE_URL + "GSE40279_average_beta.txt.gz",
                 "d980065c653381ce29b5bfa54f462a1dd99853bb3796cccd0ec15d78043eb953",
                 max_bytes=MATRIX_BYTES)
KEY = GeoFile(BASE_URL + "GSE40279_sample_key.txt.gz",
              "dbb7c510da6a90bd2ff203564288b3b3400298c9d50217377a347f2ef2a9e7f9")
METADATA = GeoFile(GEO + "&targ=gsm&form=text&view=brief",
                   "b007e0448000fb177c405cceb704dd4165a996ee9408dd22bdfdb0877b910fa8",
                   local_name="GSE40279_sample_metadata.soft")
PLATFORM_BYTES = 65 * 1024**2
PLATFORM = GeoFile("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/soft/GPL13534_family.soft.gz",
                   "c8ca5ad3c0ac18fcf2870daa63c127bcd5f6fea128f292d769ce335b3027822e",
                   max_bytes=PLATFORM_BYTES, local_name="GPL13534_platform_prefix.soft.gz")
PLATFORM_SOURCE_BYTES = 230_233_882_624
PROBE_COUNT = 473_034
AUTOSOMES = {str(number) for number in range(1, 23)}


def acquire_study(cache=Path("data/cache"), *, offline=False):
    root = cache / ACCESSION
    meta = acquire(METADATA, root, offline=offline)
    key = acquire(KEY, root, offline=offline)
    platform = acquire_large(PLATFORM, root, PLATFORM_BYTES, source_size=PLATFORM_SOURCE_BYTES,
                             offline=offline)
    matrix = acquire_large(MATRIX, root, MATRIX_BYTES, offline=offline)
    return tuple(item[0] for item in (matrix, meta, key, platform)), [
        item[1] for item in (matrix, meta, key, platform)
    ]


def parse_phenotypes(soft_path, key_path, *, expected_samples=656):
    records = []
    with soft_path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            label, _, value = line.strip().partition(" = ")
            if label == "^SAMPLE":
                records.append({"sample_id": value})
            elif records and label in {"!Sample_title", "!Sample_platform_id", "!Sample_characteristics_ch1"}:
                record = records[-1]
                if label == "!Sample_characteristics_ch1":
                    field, _, value = value.partition(": ")
                else:
                    field = label
                if field in record:
                    raise DataValidationError(f"Duplicate GEO sample field: {field}.")
                record[field] = value
    meta = pd.DataFrame(records).rename(columns={"age (y)": "age", "gender": "sex",
                                                 "!Sample_platform_id": "platform"})
    required = {"sample_id", "!Sample_title", "age", "sex", "source", "plate", "ethnicity", "tissue", "platform"}
    if len(meta) != expected_samples or not required <= set(meta.columns):
        raise DataValidationError("GSE40279 sample metadata size or schema is unexpected.")
    titles = meta["!Sample_title"].str.extract(r"^age (\d+)y (\d+)$")
    meta["age"] = pd.to_numeric(meta["age"], errors="coerce")
    meta["matrix_id"] = "X" + titles[1]
    with gzip.open(key_path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if (len(rows) != len(meta) or any(len(row) != 3 or row[0] != str(i) for i, row in enumerate(rows, 1))
            or any(not re.fullmatch(r"\d+", row[1]) or not re.fullmatch(r"\d+_R\d+C\d+", row[2])
                   for row in rows)):
        raise DataValidationError("Invalid GSE40279 sample key.")
    key = pd.DataFrame({"matrix_id": ["X" + row[1] for row in rows], "array_position": [row[2] for row in rows]})
    if (not meta.sample_id.str.fullmatch(r"GSM\d+").all() or meta.sample_id.duplicated().any()
            or meta.matrix_id.duplicated().any() or key.matrix_id.duplicated().any()
            or key.array_position.duplicated().any() or set(key.matrix_id) != set(meta.matrix_id)
            or not meta.age.between(0, 120).all()
            or not meta.age.eq(pd.to_numeric(titles[0], errors="coerce")).all()
            or not meta.sex.isin(["F", "M"]).all() or not meta.plate.str.fullmatch(r"\d+").all()
            or not meta.source.fillna("").ne("").all() or not meta.ethnicity.fillna("").ne("").all()
            or not meta.tissue.eq("whole blood").all() or not meta.platform.eq("GPL13534").all()):
        raise DataValidationError("Sample ID, sample key, age, tissue, platform or covariates do not match.")
    return meta.drop(columns="!Sample_title").merge(key, on="matrix_id", validate="one_to_one")


def parse_platform(path):
    fields = {
        "ID": "cpg_id", "Genome_Build": "build", "CHR": "chromosome", "MAPINFO": "position",
        "UCSC_RefGene_Name": "gene_symbols", "Relation_to_UCSC_CpG_Island": "island_relation",
        "Probe_SNPs_10": "nearby_snp",
    }
    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        if not any(line.strip() == "!platform_table_begin" for line in handle):
            raise DataValidationError("Missing GPL13534 platform table marker.")
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader)
        if not fields.keys() <= set(header):
            raise DataValidationError("Missing required GPL13534 annotation fields.")
        indexes = [header.index(field) for field in fields]
        for row in reader:
            if row == ["!platform_table_end"]:
                break
            if row and re.fullmatch(r"cg\d{8}", row[0]):
                if max(indexes) >= len(row):
                    raise DataValidationError("Malformed GPL13534 annotation row.")
                rows.append(tuple(row[i] for i in indexes))
        else:
            raise DataValidationError("Incomplete GPL13534 platform table.")
    anno = pd.DataFrame.from_records(rows, columns=fields.values()).set_index("cpg_id")
    if anno.empty or not anno.index.is_unique:
        raise DataValidationError("Missing or duplicate CpG annotations in GPL13534.")
    anno["position"] = pd.to_numeric(anno.position, errors="coerce").astype("Int64")
    anno["nearby_snp"] = anno.nearby_snp.ne("")
    for field in ("chromosome", "build", "island_relation"):
        anno[field] = anno[field].astype("category")
    return anno


def beta_chunks(path, metadata, *, chunk_rows=512):
    """Yield CpG IDs and numeric float64 blocks in explicit metadata order."""
    if not 1 <= chunk_rows <= 2048:
        raise ValueError("chunk_rows must be between 1 and 2048.")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        header = stream.readline().rstrip("\r\n").split("\t")
        identifiers = header[1:]
        if (header[0] != "ID_REF" or len(identifiers) != len(metadata) or len(identifiers) > 656
                or len(set(identifiers)) != len(identifiers) or set(identifiers) != set(metadata.matrix_id)):
            raise DataValidationError("Beta matrix samples do not match the GEO metadata and sample key.")
        positions = {identifier: index for index, identifier in enumerate(identifiers)}
        order = [positions[identifier] for identifier in metadata.matrix_id]
        reordered = order != list(range(len(order)))
        missing_tokens = {"", "NA", "N/A", "NULL", "null", "NaN", "nan", "#N/A", "<NA>"}
        ids = []
        block = np.empty((chunk_rows, len(order)), dtype=np.float64)
        for lineno, line in enumerate(stream, 2):
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(order) + 1 or not values[0]:
                raise DataValidationError(f"Invalid beta matrix row width or probe identifier at line {lineno}.")
            fields = values[1:]
            try:
                numeric = np.asarray(fields, dtype=np.float64)
            except ValueError:
                try:
                    numeric = np.asarray([np.nan if value in missing_tokens else value
                                          for value in fields], dtype=np.float64)
                except ValueError as exc:
                    raise DataValidationError(f"Invalid numeric beta matrix at line {lineno}.") from exc
            block[len(ids)] = numeric[order] if reordered else numeric
            ids.append(values[0])
            if len(ids) == chunk_rows:
                yield pd.Series(ids, dtype=str), block
                ids = []
                block = np.empty_like(block)
        if ids:
            yield pd.Series(ids, dtype=str), block[:len(ids)]
