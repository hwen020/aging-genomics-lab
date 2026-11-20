"""Browserless app smoke test; runs on CI when Streamlit is installed."""
from importlib.util import find_spec
from pathlib import Path

import pytest


@pytest.mark.skipif(find_spec("streamlit") is None,
                    reason="Streamlit unavailable in the offline artifact-building environment")
def test_dashboard_renders_in_streamlit_apptest():
    from streamlit.testing.v1 import AppTest

    root = Path(__file__).resolve().parents[1]
    app = AppTest.from_file(str(root / "app.py"), default_timeout=45)
    app.run()
    assert not app.exception, [x.message for x in app.exception]
