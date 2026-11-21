"""Compare Zahn Table S2 with the existing, fixed GSE5086 age analysis."""

from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from . import __version__
from .data import DataValidationError
from .ingest.geo import GeoFile, acquire
from .public_studies import ACCESSION, PAPER

TABLE = GeoFile(
    "https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=10.1371/journal.pgen.0020115.st002",
    "fb607fa3c9940f60ea94c436240888ece86fab2048f757ae9ad69219fc1c9e45",
    max_bytes=1024 * 1024, local_name="Zahn_Table_S2.xls",
)
FDR = 0.05


def read_table(path):
    raw = pd.read_excel(path, header=None, engine="xlrd")
    genes = raw.loc[raw[0].astype(str).str.fullmatch(r"Hs\.\d+|\d+(?:_[asxg])?_at")].copy()
    genes.columns = ["published_id", "published_symbol", "published_name",
                     "published_log2_slope_per_year", "published_table_probability"]
    for column in ("published_log2_slope_per_year", "published_table_probability"):
        genes[column] = pd.to_numeric(genes[column], errors="coerce")
    probability = genes.published_table_probability
    slope = genes.published_log2_slope_per_year
    tails = ((slope < 0) & probability.between(0, 0.001, inclusive="neither")
             | (slope > 0) & probability.between(0.999, 1, inclusive="neither"))
    if (len(genes) != 250 or genes.published_id.duplicated().any()
            or not tails.all()):
        raise DataValidationError("Unexpected Zahn Table S2 entries or significance values.")
    return genes.drop(columns="published_name")


def compare_symbols(published, ages):
    """Screen exact, unambiguous GPL570 symbols without choosing a best probe."""
    columns = ["gene_id", "gene_symbol", "log2_change_per_decade", "fdr"]
    mapped = ages.loc[ages.gene_symbol.notna(), columns].copy()
    mapped["symbol_key"] = mapped.gene_symbol.str.strip().str.upper()
    mapped = mapped.loc[~mapped.symbol_key.str.contains("/", regex=False)
                        & ~mapped.symbol_key.isin(["", "---", "NA"])]
    published = published.copy()
    published["symbol_key"] = published.published_symbol.fillna("").str.strip().str.upper()
    direct = published.published_id.str.endswith("_at")
    probes = published.loc[direct].merge(ages[columns], left_on="published_id", right_on="gene_id", how="left")
    symbols = published.loc[~direct].merge(mapped, on="symbol_key", how="left")
    pairs = pd.concat([probes.assign(match_rule="GPL570 probe ID"),
                       symbols.assign(match_rule="Exact gene symbol")], ignore_index=True)
    pairs["sign_agreement"] = (
        pairs.published_log2_slope_per_year.mul(pairs.log2_change_per_decade).gt(0)
        .where(pairs.gene_id.notna())
    )
    return pairs.drop(columns="symbol_key")


def compare_study(study, age, cache=Path("data/cache"), *, offline=False):
    path, table_record = acquire(TABLE, cache / ACCESSION, offline=offline)
    genes = read_table(path)
    pairs = compare_symbols(genes, age)
    matched = pairs.dropna(subset=["gene_id"])
    local_hits = int((age.fdr < FDR).sum())
    counts = {
        "published_rows": len(genes), "matched_published_rows": int(matched.published_id.nunique()),
        "unique_matched_probes": int(matched.gene_id.nunique()),
        "matched_pairs": len(matched), "same_direction_pairs": int(matched.sign_agreement.sum()),
        "local_age_fdr_hits": local_hits, "linked_age_fdr_hits": int((matched.fdr < FDR).sum()),
    }
    report = f"""# Published Study Comparison · GSE5086

## Zahn et al. (2006)

The paper reported **250 of 31,948 genes** at nominal age p < 0.001. Age permutation
with sex and anatomy held with each sample yielded about 32 hits by chance on average, an
estimated false-positive fraction of about 13% among the 250. The model used all 81 muscle
samples and least squares with continuous age, sex and a two-level anatomy covariate.
dChip PM-only probe-set intensities were averaged by gene on the linear scale before log2.
In the Table S2 workbook, values in its p-value column fall below 0.001 for negative
slopes and above 0.999 for positive slopes. The export preserves these as published
table probabilities and does not treat them as two-sided p-values.

## Current P1 age analysis

The cohort includes {len(study[0].metadata)} abdominal samples and {age.p_value.count():,} tested
GPL570 probes. The test uses continuous age, sex, HC3 errors and BH across valid probe
p-values. There were **{local_hits} probe-level FDR < 0.05 hits**. P1 uses the deposited
dChip scale, filters nonpositive or nonfinite probe rows and applies log2 at probe level.
Its separate age-group contrast does not enter this comparison.

## Table S2 identifier and direction audit

{counts['matched_published_rows']} of 250 published rows matched at least one retained
probe, using the original GPL570 probe ID when available and otherwise using an
exact case-insensitive symbol with a single GPL570 symbol. The links involve
{counts['unique_matched_probes']} distinct GPL570 probes and form
{counts['matched_pairs']} published-record–probe pairs. In {counts['same_direction_pairs']} pairs
the estimated age-effect direction agreed. {counts['linked_age_fdr_hits']} pairs include
a P1 probe with FDR < 0.05. The direction count is descriptive because both analyses
reuse samples and pairs may share published records or probes.

## Remaining uncertainty

The difference combines cohort, mapping, aggregation, standard errors and the significance
rule. Their individual contributions have not been isolated. Historical UniGene identifiers
and current symbols may differ, and an unmatched row does not imply biological disagreement.
The nominal paper threshold and its permutation estimate are distinct from a per-probe
BH threshold. The primary P1 models and decision thresholds were used as originally specified.

[Paper]({PAPER}) · [Published Table S2](https://doi.org/10.1371/journal.pgen.0020115.st002)
"""
    provenance = {"version": __version__, "run_utc": datetime.now(timezone.utc).isoformat(),
                  "xlrd": version("xlrd"), "paper_url": PAPER,
                  "source_files": [*study[4], table_record], "local_fdr": FDR, "counts": counts}
    return report, pairs, provenance
