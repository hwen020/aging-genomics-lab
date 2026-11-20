import json
from pathlib import Path

import pandas as pd

from aginggenomics.reporting import export_analysis, sha256
from aginggenomics.simulation import write_demo


def test_export_has_provenance_and_results(tmp_path: Path):
    data = tmp_path / "demo"
    out = tmp_path / "results"
    write_demo(data)
    summary = export_analysis(data / "expression.csv", data / "metadata.csv", out,
                              simulated=True, source="test-synthetic")
    assert summary["qc"]["n_samples"] == 60
    assert summary["expression_sha256"] == sha256(data / "expression.csv")
    for name in ("differential_expression.csv", "age_association.csv",
                 "sample_pca.csv", "report.md", "provenance.json"):
        assert (out / name).exists()
    assert len(pd.read_csv(out / "differential_expression.csv")) == 320
    assert "SYNTHETIC" in (out / "report.md").read_text(encoding="utf-8")
    assert json.loads((out / "provenance.json").read_text())["input_type"] == "synthetic"


def test_demo_manifest_truth(tmp_path: Path):
    write_demo(tmp_path)
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["human_subjects"] is False
    assert set(pd.read_csv(tmp_path / "simulation_truth.csv").injected_pattern) == {
        "group_shift", "age_trend", "background"
    }
