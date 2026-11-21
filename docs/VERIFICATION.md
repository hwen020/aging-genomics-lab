# Local verification · Windows · 2026-11-20

## v0.3.0 acceptance

The extracted v0.2.0 workspace had no recoverable Git metadata. Local baseline commit `2219392` includes the reviewed source, documentation and bundled synthetic demo. The pre-change test run passed 39/39. No remote is configured.

```powershell
.venv\Scripts\python.exe -m pip install -e ".[dev]"
.venv\Scripts\python.exe -m pip check
.venv\Scripts\python.exe -m pytest
.venv\Scripts\python.exe -m compileall -q src app.py scripts/fetch_public_dataset.py
.venv\Scripts\python.exe -m aginggenomics.cli public --offline --out output/GSE5086-P2-check
.venv\Scripts\python.exe -m aginggenomics.cli compare-published --offline
.venv\Scripts\python.exe -m aginggenomics.cli ewas --offline
.venv\Scripts\python.exe -m aginggenomics.cli ewas --offline --out output/GSE40279-repeat
scripts\start_windows.bat
curl.exe -fsS http://127.0.0.1:8501/_stcore/health
```

- The seven P1 CSVs are byte-identical to the v0.2.0 run. The independent comparison report links 113 published rows to 219 probe pairs with 196 direction agreements.
- The GSE40279 matrix was downloaded in full from GEO. The verified cache contains a SHA-256 pinned annotation prefix from the official platform table. A separate live NCBI smoke test resumed a 101-byte bounded artifact from 40 local bytes using an exact 61-byte range.
- Both offline EWAS runs analyzed 656 samples and generated 424,870 tests and 58,849 BH-FDR hits. All seven output CSVs are byte-identical across runs. On this Windows host, the measured repeat took 134.3 seconds with a peak working set of 601.1 MiB, using Python 3.11.7 and two OpenBLAS threads.
- The final suite passed **51/51** locally, with 12 new tests. Dependency and compilation checks passed. Chromium verified the Synthetic Demo, both Public Study workflows, loading and rerunning the full EWAS from verified local cache, the PCA/volcano/CpG plots, CpG selection and search, CSV and ZIP downloads, and the published-study comparison. [EWAS overview](screenshots/ewas-overview.png) · [EWAS analysis](screenshots/ewas-analysis.png)
- Git ignore verification covers raw/cache artifacts, exports, CSV uploads and representative credentials. GitHub Actions is configured to compile, test with small fixtures and run the synthetic CLI on Python 3.10 and 3.12. No hosted workflow was run for this checkout.

## Archived P1 acceptance

## Workspace baseline

The supplied folder was an extracted v0.1.0 workspace. `git status --short --branch` and `git log -1 --oneline` found no Git metadata. The baseline suite passed **26/26** under Python 3.11.7, and the baseline CLI exported 60 synthetic samples.

| Pre-change file | SHA-256 |
|---|---|
| `app.py` | `168808d80ea7bf6ebc6a1484fcf604772bcb56fc7aebad9d0eebb438fc167cc4` |
| `README.md` | `968a4c9788e9b88097268191f6ef74cad86e431241f25b2882a87dd5e3eb18ae` |
| `pyproject.toml` | `4a822e6ee2f3661f76813f7c5ec193dadac75b5f98042809a9c8db214af7afc8` |

## P1 acceptance commands

```powershell
.venv\Scripts\python.exe -m pip install -e ".[dev]"
.venv\Scripts\python.exe -m pip check
.venv\Scripts\python.exe -m pytest
.venv\Scripts\python.exe -m compileall -q src app.py scripts/fetch_public_dataset.py
.venv\Scripts\python.exe -m aginggenomics.cli analyze --expression data/demo/expression.csv --metadata data/demo/metadata.csv --out output/demo-analysis --simulated --source SYNTHETIC_AGING_DEMO_V1
.venv\Scripts\python.exe -m aginggenomics.cli public --cache data/cache/acceptance --out output/GSE5086
.venv\Scripts\python.exe -m aginggenomics.cli public --offline --out output/GSE5086
.venv\Scripts\python.exe scripts/fetch_public_dataset.py --offline
scripts\start_windows.bat
curl.exe -fsS http://127.0.0.1:8501/_stcore/health
```

The final test suite passed **39/39** on Python 3.11.7 with the declared dependencies. Compilation and `pip check` passed. The fresh-cache public run verified both official downloads and generated a report for 62 samples and 54,615 tested probes. Comparing two offline exports from independent verified caches gave identical SHA-256 hashes for all seven CSVs. The synthetic CLI exported 60 samples and 320 features using bundled CSVs.

The Windows batch launcher started Streamlit, and its health endpoint returned `ok`. A headless Chromium session at 1440 × 960 exercised Synthetic Demo, all four public tabs, cached-data verification, five interactive charts, probe search and both CSV and ZIP downloads. The downloaded ZIP contained the report, provenance and chart artifacts. [Overview](screenshots/public-overview.png) · [Analysis](screenshots/public-analysis.png) · [Results](screenshots/public-results.png)

The GitHub Actions definition runs fixture-based tests and a synthetic CLI smoke test. A hosted GitHub Actions run requires a Git checkout and remote.

## v0.1.0 packaging verification

- Python syntax compilation: PASS (`python -m compileall -q src app.py`).
- `pytest` core and reporting suite: PASS (25 tests).
- CLI demo regeneration: PASS (generated a deterministic synthetic dataset).
- CLI report smoke: PASS (60 simulated samples / 320 features; all five outputs created).
- Streamlit integration test: **not executed in the artifact-building environment** because `streamlit` was not installed and package download could not reach PyPI. The app syntax compiles, but **the actual browser rendering remains to be checked** after the user installs dependencies in VS Code.

The `tests/test_app.py` test uses Streamlit's `AppTest` and is designed to run automatically once dependencies are installed (including on GitHub Actions). Do not report the browser UI as fully verified until this test and a browser walkthrough pass.
