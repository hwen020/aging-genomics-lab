"""Create the aggregate-only Portfolio Demo snapshot from local release reports."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from aginggenomics import __version__
from aginggenomics.reporting import sha256


def read_study(root, accession, filename, feature, effect, hits, formula):
    folder = root / accession
    provenance = json.loads((folder / "provenance.json").read_text(encoding="utf-8"))
    if provenance["accession"] != accession or provenance["project_version"] != __version__:
        raise ValueError(f"{accession} requires a verified export from v{__version__}.")
    values, examples, hit_count = [], [], 0
    result = folder / filename
    with result.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if not row["p_value"]:
                continue
            p = float(row["p_value"])
            fdr = float(row["fdr"])
            values.append(p)
            hit_count += fdr < 0.05
            if len(examples) < 8:
                examples.append({"feature_id": row[feature], "effect_per_decade": float(row[effect]),
                                 "p_value": p, "fdr": fdr})
    counts, _ = np.histogram(values, bins=np.linspace(0, 1, 21))
    qc = provenance["qc"]
    if hit_count != hits(qc):
        raise ValueError(f"{accession} discovery count differs from provenance.")
    return {
        "accession": accession, "source": provenance["geo_url"], "paper": provenance["paper_url"],
        "generated_utc": provenance["run_utc"], "analysis_version": __version__,
        "samples": qc.get("n_samples", qc.get("samples")), "tests": len(values),
        "fdr_hits": hit_count, "fdr_threshold": 0.05,
        "formula": provenance[formula], "limitations": provenance["limitations"],
        "p_histogram": counts.tolist(), "examples": examples, "result_sha256": sha256(result),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", type=Path, default=Path("output/release"))
    parser.add_argument("--out", type=Path, default=Path("src/aginggenomics/portfolio_studies.json"))
    args = parser.parse_args()
    root = args.reports
    comparison = json.loads((root / "GSE5086-comparison/provenance.json").read_text(encoding="utf-8"))
    if comparison["version"] != __version__:
        raise ValueError("Published comparison must come from this release.")
    expression = read_study(root, "GSE5086", "age_association.csv", "gene_id",
                            "log2_change_per_decade", lambda qc: comparison["counts"]["local_age_fdr_hits"],
                            "age_formula")
    expression["published_comparison"] = comparison["counts"]
    methylation = read_study(root, "GSE40279", "ewas_results.csv", "cpg_id",
                             "m_change_per_decade", lambda qc: qc["fdr_hits"], "parameters")
    methylation["formula"] = methylation["formula"]["model_formula"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"GSE5086": expression, "GSE40279": methylation},
                                   indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
