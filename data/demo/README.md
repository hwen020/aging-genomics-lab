# Synthetic demo dataset

All 60 samples and 320 feature IDs are computer-generated from a fixed seed.
These are not public GEO data or actual human samples. See `manifest.json`.

- 30 `Younger` ages 22–40, 30 `Older` ages 62–86.
- 24 genes with injected group shifts, 12 with injected age trends.
- Some associations are confounded because age ranges define the two groups.
- `simulation_truth.csv` documents construction and supports test assertions.
- Expression values are simulated on a continuous log2-like scale and are
  **not** a realistic model of count sampling, sequencing noise or biology.

Regenerate with:

```bash
python -m aginggenomics.cli demo --out data/demo
```
