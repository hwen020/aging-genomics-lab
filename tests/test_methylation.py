import gzip
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from aginggenomics.cli import main
from aginggenomics.data import DataValidationError
from aginggenomics.ingest import resumable
from aginggenomics.ingest.geo import GeoFile
from aginggenomics.methylation_analysis import age_design, analyze_matrix
from aginggenomics.methylation_reporting import export_ewas
from aginggenomics.methylation_study import beta_chunks, parse_phenotypes, parse_platform


@pytest.fixture
def methylation_fixture(tmp_path):
    sample_ids = [f"GSM{i}" for i in range(1, 13)]
    age = [25, 31, 43, 57, 62, 80, 27, 33, 45, 58, 69, 83]
    sources = []
    key_rows = []
    for i, (sample, years) in enumerate(zip(sample_ids, age), 1):
        sources += [f"^SAMPLE = {sample}", f"!Sample_title = age {years}y {1000+i}",
                    "!Sample_platform_id = GPL13534"]
        fields = {"age (y)": years, "gender": "MF"[i % 2], "source": "Site",
                  "plate": (i - 1) // 6 + 1, "ethnicity": "Recorded", "tissue": "whole blood"}
        sources += [f"!Sample_characteristics_ch1 = {k}: {v}" for k, v in fields.items()]
        key_rows.append(f"{i}\t{1000+i}\t1000{i}_R01C01")
    soft = tmp_path / "sample.soft"
    soft.write_text("\n".join(sources), encoding="utf-8")
    key = tmp_path / "key.gz"
    with gzip.open(key, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(key_rows))
    metadata = parse_phenotypes(soft, key, expected_samples=12)
    ids = [f"cg{i*64:08d}" for i in range(1, 9)]
    mapping = tmp_path / "platform.gz"
    with gzip.open(mapping, "wt", encoding="utf-8") as handle:
        handle.write("!platform_table_begin\nID\tGenome_Build\tCHR\tMAPINFO\tUCSC_RefGene_Name"
                     "\tRelation_to_UCSC_CpG_Island\tProbe_SNPs_10\n")
        for i, probe in enumerate(ids[:-1]):
            chrom = "X" if i == 6 else "1"
            snp = "rs1" if i == 5 else ""
            handle.write(f"{probe}\t37\t{chrom}\t{1000+i}\tGENE{i}\tIsland\t{snp}\n")
        handle.write("!platform_table_end\n")
    annotation = parse_platform(mapping)
    ages = metadata.age.to_numpy()
    rng = np.random.default_rng(12)
    values = 0.2 + 0.006 * (ages - min(ages)) + rng.normal(0, 0.02, (8, 12))
    values[2, 1] = np.nan
    values[3, 1] = 1.2
    values[4] = 0.3
    matrix = tmp_path / "beta.gz"
    matrix_frame = pd.DataFrame(values, index=pd.Index(ids, name="ID_REF"), columns=metadata.matrix_id)
    matrix_frame.iloc[:, ::-1].to_csv(matrix, sep="\t", compression="gzip")
    return matrix, metadata, annotation, matrix_frame, soft, key


def test_blockwise_ewas_qc_matching_and_report(methylation_fixture, tmp_path):
    matrix, metadata, annotation, frame, _, _ = methylation_fixture
    identifiers, values = next(beta_chunks(matrix, metadata, chunk_rows=2))
    np.testing.assert_allclose(values, frame.iloc[:2].to_numpy())
    small = analyze_matrix(matrix, metadata, annotation, chunk_rows=2, expected_rows=8)
    large = analyze_matrix(matrix, metadata, annotation, chunk_rows=8, expected_rows=8)
    pd.testing.assert_frame_equal(small.results, large.results, rtol=1e-10)
    assert len(identifiers) == 2
    assert small.results.m_change_per_decade.gt(0).all()
    assert small.results.fdr.lt(0.05).all()
    assert small.pca.pc1.notna().all()
    assert small.sample_qc.missing_values.sum() == 1
    assert small.sample_qc.out_of_range_values.sum() == 1
    assert set(small.probe_qc.exclusion_reason) == {
        "", "Missing beta", "Out-of-range or infinite beta", "Low beta variability",
        "Nearby annotated SNP", "Non-autosomal probe", "Missing hg19 annotation",
    }
    out = tmp_path / "report"
    export_ewas(small, [{"sha256": "test-only-fixture"}], out)
    record = json.loads((out / "provenance.json").read_text(encoding="utf-8"))
    assert record["covariates"] == ["sex", "plate"]
    assert record["qc"]["tested_probes"] == 2
    assert (out / "chart_pca.html").exists()
    assert len(pd.read_csv(out / "significant_sites.csv")) == 2


def test_matrix_id_mismatch_and_memory_bound(methylation_fixture):
    matrix, metadata, _, _, _, _ = methylation_fixture
    with pytest.raises(DataValidationError, match="samples do not match"):
        next(beta_chunks(matrix, metadata.assign(matrix_id="unrelated"), chunk_rows=1))
    with pytest.raises(ValueError, match="2048"):
        next(beta_chunks(matrix, metadata, chunk_rows=2049))


def test_truncated_gzip_fails_instead_of_returning_completed_run(methylation_fixture):
    matrix, metadata, _, _, _, _ = methylation_fixture
    matrix.write_bytes(matrix.read_bytes()[:-8])
    with pytest.raises(EOFError):
        list(beta_chunks(matrix, metadata, chunk_rows=2))


@pytest.mark.parametrize(("replacement", "message"), [("0.5", "row width"), ("\tbogus", "Invalid numeric")])
def test_beta_row_parser_fails_on_bad_width_and_tokens(methylation_fixture, replacement, message):
    matrix, metadata, _, _, _, _ = methylation_fixture
    with gzip.open(matrix, "rt") as handle:
        lines = handle.readlines()
    lines[1] = lines[1].rstrip("\r\n").rsplit("\t", 1)[0] + replacement + "\n"
    with gzip.open(matrix, "wt") as handle:
        handle.writelines(lines)
    with pytest.raises(DataValidationError, match=message):
        list(beta_chunks(matrix, metadata, chunk_rows=2))


def test_rejects_excess_matrix_columns_before_allocating_a_block(tmp_path):
    matrix = tmp_path / "oversized.gz"
    identifiers = [f"X{i}" for i in range(657)]
    with gzip.open(matrix, "wt") as handle:
        handle.write("\t".join(["ID_REF", *identifiers]) + "\n")
    with pytest.raises(DataValidationError, match="samples do not match"):
        next(beta_chunks(matrix, pd.DataFrame({"matrix_id": identifiers})))


def test_duplicate_probe_across_blocks_is_rejected(methylation_fixture):
    matrix, metadata, annotation, frame, _, _ = methylation_fixture
    pd.concat([frame, frame.iloc[:1]]).to_csv(matrix, sep="\t", compression="gzip")
    with pytest.raises(DataValidationError, match="duplicate probe"):
        analyze_matrix(matrix, metadata, annotation, chunk_rows=3, expected_rows=9)


def test_sample_key_mismatch_and_singular_model(methylation_fixture):
    _, metadata, _, _, soft, key = methylation_fixture
    with gzip.open(key, "rt") as handle:
        text = handle.read().replace("\t1001\t", "\t9999\t")
    with gzip.open(key, "wt") as handle:
        handle.write(text)
    with pytest.raises(DataValidationError, match="do not match"):
        parse_phenotypes(soft, key, expected_samples=12)
    with pytest.raises(DataValidationError, match="Rank-deficient"):
        age_design(metadata.assign(sex="F"))


class Response(io.BytesIO):
    status = 206

    def __init__(self, payload, header):
        super().__init__(payload)
        self.headers = {"Content-Range": header}


def test_resume_after_short_response_and_offline_corruption(monkeypatch, tmp_path):
    contents = b"abcdef"
    artifact = GeoFile("https://example.invalid/test.bin", hashlib.sha256(contents).hexdigest())
    responses = iter([Response(b"abc", "bytes 0-5/6"), Response(b"def", "bytes 3-5/6")])
    requested = []

    def get(request, timeout):
        requested.append(request.headers["Range"])
        return next(responses)

    monkeypatch.setattr(resumable, "urlopen", get)
    monkeypatch.setattr(resumable.time, "sleep", lambda _seconds: None)
    path, _ = resumable.acquire_large(artifact, tmp_path, 6)
    assert requested == ["bytes=0-5", "bytes=3-5"]
    assert path.read_bytes() == contents
    path.write_bytes(b"ABCDEF")
    with pytest.raises(OSError, match="No verified cache"):
        resumable.acquire_large(artifact, tmp_path, 6, offline=True)


def test_reject_wrong_range_and_budget(monkeypatch, tmp_path):
    artifact = GeoFile("https://example.invalid/test.bin", "0" * 64, max_bytes=6)
    monkeypatch.setattr(resumable, "urlopen", lambda *a, **k: Response(b"", "bytes 0-400/401"))
    monkeypatch.setattr(resumable.time, "sleep", lambda _seconds: None)
    with pytest.raises(OSError, match="Content-Range"):
        resumable.acquire_large(artifact, tmp_path, 6, source_size=401)
    assert not (tmp_path / artifact.filename).exists()
    with pytest.raises(ValueError, match="transfer limit"):
        resumable.acquire_large(artifact, tmp_path, 7)


def test_ewas_cli_validates_chunk_limit(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["aginggenomics", "ewas", "--chunk-rows", "4096"])
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 2
    assert "2048" in capsys.readouterr().err


def test_dashboard_selects_methylation_without_downloading():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_file(Path(__file__).resolve().parents[1] / "app.py", default_timeout=45).run()
    app.sidebar.radio[0].set_value("Public Study").run()
    app.selectbox[0].set_value("GSE40279").run()
    assert not app.exception
    assert "Methylation Analysis" in [item.value for item in app.header]
