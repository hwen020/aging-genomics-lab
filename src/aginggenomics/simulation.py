"""Seeded, labeled synthetic test data; never representative of real subjects."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def make_demo(seed: int = 20261119, n_genes: int = 320,
              per_group: int = 30) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    n = 2 * per_group
    groups = np.array(["Younger"] * per_group + ["Older"] * per_group)
    ages = np.r_[rng.integers(22, 41, size=per_group), rng.integers(62, 87, size=per_group)]
    samples = [f"SIM_{i + 1:03d}" for i in range(n)]
    genes = [f"SIM_GENE_{i+1:04d}" for i in range(n_genes)]
    batch = np.array(["B1" if i % 2 else "B2" for i in range(n)])
    sex = np.array(["F" if i % 2 else "M" for i in range(n)])
    # Synthetic Gaussian log2-scale expression (not integer RNA-seq counts).
    baseline = rng.uniform(5, 10.5, size=(n_genes, 1))
    expression = baseline + rng.normal(0, 0.72, size=(n_genes, n))
    n_diff = min(24, n_genes // 5)
    n_age = min(12, n_genes // 8)
    delta = np.r_[rng.uniform(1.0, 1.8, size=n_diff // 2),
                  -rng.uniform(1.0, 1.8, size=n_diff - n_diff // 2)]
    expression[:n_diff, per_group:] += delta[:, None]
    # Independent demonstration signal; group/age remain confounded by design.
    age_effects = rng.choice([-1, 1], size=n_age) * rng.uniform(0.18, 0.36, size=n_age)
    expression[n_diff:n_diff + n_age] += age_effects[:, None] * ((ages - ages.mean()) / 10)
    expr = pd.DataFrame(expression, columns=samples)
    expr.insert(0, "gene_id", genes)
    meta = pd.DataFrame({"sample_id": samples, "group": groups, "age": ages,
                         "batch": batch, "sex": sex})
    truth = pd.DataFrame({
        "gene_id": genes,
        "injected_pattern": (["group_shift"] * n_diff + ["age_trend"] * n_age +
                             ["background"] * (n_genes - n_diff - n_age)),
    })
    return expr, meta, truth


def write_demo(output_dir: str | Path) -> None:
    import json
    p = Path(output_dir)
    p.mkdir(parents=True, exist_ok=True)
    expr, meta, truth = make_demo()
    expr.to_csv(p / "expression.csv", index=False, float_format="%.8f")
    meta.to_csv(p / "metadata.csv", index=False)
    truth.to_csv(p / "simulation_truth.csv", index=False)
    (p / "manifest.json").write_text(json.dumps({
        "dataset_id": "SYNTHETIC_AGING_DEMO_V1",
        "data_origin": "synthetic",
        "human_subjects": False,
        "method": "numpy seeded Gaussian log2 scale, group shifts and age trends",
        "seed": 20261119,
        "meaning": "Software validation and presentation ONLY; no real biological claims",
        "expression_scale": "normalized_log2_like_simulated_values",
        "license": "CC0-1.0",
    }, indent=2) + "\n", encoding="utf-8")
