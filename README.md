# DESI Y1 Analysis for Dark Geometry

This directory contains scripts to analyze DESI Year 1 data with Dark Geometry.

## Files

| File | Description |
|------|-------------|
| `desi_y1_analysis.py` | Main comparison: DG vs ΛCDM vs w₀wₐCDM |
| `w_evolution.py` | Dark energy equation of state w(z) |
| `desi_likelihood.py` | Likelihood for MCMC analysis |

## Quick Start

```bash
# Run main analysis
python desi_y1_analysis.py

# Analyze w(z) evolution
python w_evolution.py

# Quick likelihood comparison
python desi_likelihood.py
```

## DESI Y1 Data

### BAO Measurements (arXiv:2404.03000)

| Tracer | z_eff | D_M/r_d | D_H/r_d |
|--------|-------|---------|---------|
| BGS | 0.295 | — | D_V/r_d = 7.93 |
| LRG1 | 0.510 | 13.62 | 20.98 |
| LRG2 | 0.706 | 16.85 | 20.08 |
| LRG3 | 0.930 | 21.71 | 17.88 |
| ELG | 1.317 | 27.79 | 13.82 |
| QSO | 1.491 | 30.69 | 13.26 |
| Lyα | 2.330 | 39.71 | 8.52 |

### RSD Measurements

| z_eff | f·σ₈ |
|-------|------|
| 0.295 | 0.392 ± 0.044 |
| 0.510 | 0.458 ± 0.033 |
| 0.706 | 0.449 ± 0.032 |
| 0.930 | 0.437 ± 0.035 |
| 1.317 | 0.372 ± 0.062 |

### w₀-wₐ Constraints (arXiv:2404.03002)

DESI + CMB:
- w₀ = -0.55 ± 0.21
- wₐ = -1.30 ± 0.70

## Key Results

### DG Predictions

| Parameter | DG Value | DESI |
|-----------|----------|------|
| w₀ | -0.88 | -0.55 ± 0.21 |
| wₐ | -0.13 | -1.30 ± 0.70 |
| σ₈ | 0.773 | — |

### Why DESI Matters for DG

1. **Dynamic w(z)**: DESI finds w₀ > -1, consistent with DG transition
2. **f·σ₈**: DG predicts lower values, matches trend
3. **No fine-tuning**: DG derives w(z) from first principles

## Future Forecasts

| Release | σ(w₀) | σ(f·σ₈) |
|---------|-------|---------|
| Y1 (2024) | 0.21 | 0.03-0.06 |
| Y3 (2026) | ~0.12 | ~0.02 |
| Y5 (2028) | ~0.09 | ~0.015 |

With Y5 data, DG vs ΛCDM discrimination will be **decisive** (Δχ² ~ 25-50).

## References

1. DESI 2024 III: BAO (arXiv:2404.03000)
2. DESI 2024 IV: Cosmology (arXiv:2404.03001)
3. DESI 2024 V: Full-shape (arXiv:2404.03002)
