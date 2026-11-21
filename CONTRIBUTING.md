# Contributing

Use Python 3.10 or later and install with `python -m pip install -e ".[dev]"`. Run `python -m pytest` from the repository root. Prefer small changes with fixture-based tests that run without downloading GEO data.

Document any change to sample selection, probe filtering, effect coding, model assumptions or multiple testing. Regenerate full reports locally before updating a portfolio snapshot. Keep the results and the versioned snapshot traceable to their accession, source checksums and analysis version.

Keep raw or participant-level records, credentials, machine-specific paths and full reports out of commits. Public examples should contain only aggregate counts and feature-level estimates. Use the synthetic fixture for debugging and bug reports.

Contributions to project code are accepted under the repository's [MIT license](LICENSE). Cite the project and original dataset sources. See [source terms](docs/PUBLIC_STUDIES.md).
