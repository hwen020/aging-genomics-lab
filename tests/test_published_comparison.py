import pandas as pd
import pytest

from aginggenomics.data import DataValidationError
from aginggenomics.published_comparison import compare_symbols, read_table


def test_direct_and_strict_symbol_matching_retains_all_eligible_probes():
    paper = pd.DataFrame({"published_id": ["Hs.1", "123_at", "Hs.3"],
                          "published_symbol": ["GeneA", "---", "GeneB"],
                          "published_log2_slope_per_year": [0.1, -0.2, 0.3]})
    probes = pd.DataFrame({"gene_id": ["101_at", "102_at", "123_at", "789_at"],
                           "gene_symbol": ["GENEA", "GENEA", None, "GENEB///ANOTHER"],
                           "log2_change_per_decade": [0.5, -0.2, -0.3, 1.0],
                           "fdr": [0.6, 0.7, 0.8, 0.9]})
    pairs = compare_symbols(paper, probes)
    assert len(pairs) == 4
    assert pairs.loc[pairs.published_id.eq("123_at"), "sign_agreement"].item()
    assert pairs.loc[pairs.published_id.eq("Hs.3"), "gene_id"].isna().all()
    assert pairs.loc[pairs.published_id.eq("Hs.1"), "gene_id"].tolist() == ["101_at", "102_at"]


def test_table_validates_tail_direction_without_inventing_two_sided_p(monkeypatch):
    content = pd.DataFrame([[f"Hs.{i}", "X", "Name", (-1 if i % 2 else 1),
                             (0.0005 if i % 2 else 0.9995)] for i in range(250)])
    monkeypatch.setattr(pd, "read_excel", lambda *a, **k: content)
    parsed = read_table("local-fixture.xls")
    assert parsed.published_table_probability.max() == 0.9995
    content.loc[0, 4] = 0.2
    with pytest.raises(DataValidationError, match="significance"):
        read_table("local-fixture.xls")
