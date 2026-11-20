import csv
import gzip
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from aginggenomics.cli import main
from aginggenomics.data import DataValidationError
from aginggenomics.ingest import geo
from aginggenomics.public_analysis import analyze_study, build_export, fit_probes
from aginggenomics.public_studies import OTHER_TISSUE, TISSUE, parse_annotation, parse_matrix, prepare_study


@pytest.fixture
def matrix_rows():
    samples = [f"GSM{i}" for i in range(1, 9)]
    return [
        ["!Series_geo_accession", "GSE5086"],
        ["!Sample_geo_accession", *samples],
        ["!Sample_platform_id", *(["GPL570"] * 8)],
        ["!Sample_data_processing", *(["dChip PM only"] * 8)],
        ["!Sample_source_name_ch1", *([TISSUE] * 7), OTHER_TISSUE],
        ["!Sample_characteristics_ch1", *[f"Age(yrs.): {age}" for age in [25, 28, 35, 42, 66, 70, 81, 52]]],
        ["!Sample_characteristics_ch1", *["Sex: Female", "Sex: Male"] * 4],
        ["!series_matrix_table_begin"],
        ["ID_REF", *samples[::-1]],
        ["p_keep", 345.5, 24.5, 20, 16, 12, 10, 9, 7],
        ["p_jitter", 4, 8.9, 3, 9, 1, 5, 10, 4],
        ["p_zero", 15.5, 16, 17, 18, 0, 20, 21, 22],
        ["p_missing", 15.5, 16, 17, 18, "", 20, 21, 22],
        ["p_constant", *([8] * 8)],
        ["!series_matrix_table_end"],
    ]


def write_matrix(tmp_path, rows):
    path = tmp_path / "matrix.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as stream:
        csv.writer(stream, delimiter="\t").writerows(rows)
    return path


def test_parse_align_filter_and_export(tmp_path, matrix_rows):
    expression, metadata = parse_matrix(write_matrix(tmp_path, matrix_rows))
    ds, samples, filters = prepare_study(expression, metadata)
    assert ds.expression.columns.tolist() == [f"GSM{i}" for i in range(1, 8)]
    assert ds.expression.loc["p_keep", "GSM1"] == pytest.approx(np.log2(7))
    assert samples.loc[7, "exclusion_reason"] == "Anatomical source outside rectus abdominis"
    assert filters.set_index("gene_id").exclusion_reason.to_dict() == {
        "p_keep": "", "p_jitter": "", "p_zero": "Zero intensity",
        "p_missing": "Nonfinite signal", "p_constant": "Near-constant log2 signal",
    }
    annotation = tmp_path / "annotation.gz"
    with gzip.open(annotation, "wt", encoding="utf-8") as stream:
        stream.write("^PLATFORM = GPL570\n#ID = Probe\n!platform_table_begin\n"
                     "ID\tGene symbol\tGene ID\np_keep\tA///B\t1///2\n!platform_table_end\n")
    study = (ds, samples, filters, parse_annotation(annotation), [{"sha256": "test_only"}])
    tables = analyze_study(study)
    age = tables["age_association"].set_index("gene_id")
    assert age.loc["p_keep", "log2_change_per_decade"] > 0
    assert age.loc["p_keep", "gene_symbol"] == "A///B"
    assert age.fdr.between(0, 1).all()
    files, provenance = build_export(study, tables)
    assert {"report.md", "report.html", "chart_pca.html", "probe_filters.csv"} <= files.keys()
    assert json.loads(files["provenance.json"])["covariates"] == ["sex"]
    assert provenance["qc"]["n_samples"] == 7


@pytest.mark.parametrize(("row", "col", "replacement", "error"), [
    (8, 1, "GSM999", "IDs differ"),
    (8, 1, "GSM7", "IDs differ"),
    (5, 1, "Age(yrs.):", "missing age"),
    (3, 1, "raw counts", "Counts"),
    (9, 1, "unknown", "Non-numeric"),
])
def test_invalid_series_rejected(tmp_path, matrix_rows, row, col, replacement, error):
    matrix_rows[row][col] = replacement
    with pytest.raises(DataValidationError, match=error):
        parse_matrix(write_matrix(tmp_path, matrix_rows))


def test_integer_signals_rejected_and_rank_checked(tmp_path, matrix_rows):
    expression, metadata = parse_matrix(write_matrix(tmp_path, matrix_rows))
    with pytest.raises(DataValidationError, match="continuous"):
        prepare_study(expression.round(), metadata)
    ds, _, _ = prepare_study(expression, metadata)
    with pytest.raises(DataValidationError, match="Rank-deficient"):
        fit_probes(ds, ds.metadata.sex.eq("Male"), "effect")


def test_pinned_download_retry_and_offline_corruption(tmp_path, monkeypatch):
    content = b"test-only GEO fixture"
    artifact = geo.GeoFile("https://example.invalid/test.gz", hashlib.sha256(content).hexdigest())
    responses = iter([OSError("network down"), io.BytesIO(content)])

    def open_fixture(*args, **kwargs):
        response = next(responses)
        if isinstance(response, OSError):
            raise response
        return response

    monkeypatch.setattr(geo, "urlopen", open_fixture)
    monkeypatch.setattr(geo.time, "sleep", lambda _: None)
    path, manifest = geo.acquire(artifact, tmp_path)
    assert path.read_bytes() == content
    assert geo.acquire(artifact, tmp_path, offline=True) == (path, manifest)
    path.write_bytes(b"x" * len(content))
    with pytest.raises(OSError, match="No verified cache"):
        geo.acquire(artifact, tmp_path, offline=True)


@pytest.mark.parametrize(("limit", "message"), [(1, "size"), (200, "SHA-256")])
def test_rejected_download_never_publishes_cache(tmp_path, monkeypatch, limit, message):
    artifact = geo.GeoFile("https://example.invalid/test.gz", "0" * 64, max_bytes=limit)
    monkeypatch.setattr(geo, "urlopen", lambda *a, **k: io.BytesIO(b"not the expected bytes"))
    monkeypatch.setattr(geo.time, "sleep", lambda _: None)
    with pytest.raises(OSError, match=message):
        geo.acquire(artifact, tmp_path)
    assert not list(tmp_path.iterdir())


def test_cli_blocks_declared_counts(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["aginggenomics", "analyze", "--expression", "x.csv",
                                     "--metadata", "m.csv", "--out", str(tmp_path),
                                     "--data-type", "rna-seq-counts", "--ack-normalized-log2"])
    with pytest.raises(SystemExit, match="2"):
        main()
    assert "count-aware" in capsys.readouterr().err
    assert not list(tmp_path.iterdir())


def test_dashboard_can_select_public_catalog():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_file(Path(__file__).resolve().parents[1] / "app.py", default_timeout=45).run()
    app.sidebar.radio[0].set_value("Public Study").run()
    assert not app.exception
    assert [tab.label for tab in app.tabs] == ["Overview", "Data Explorer", "Analysis", "Results & Reports"]
    app.selectbox[0].set_value("GSE164471").run()
    assert not app.exception
    assert "count-aware" in app.info[0].value
