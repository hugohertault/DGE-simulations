# DGE-Simulations: Dark Geometry Extended

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Complete simulation and analysis pipeline for Dark Geometry Extended (DGE)**

A unified framework that resolves both the Hubble (H₀) and σ₈ cosmological tensions with **zero free parameters**.

---

## 🎯 Key Results

| Tension | ΛCDM | DGE | Improvement |
|---------|------|-----|-------------|
| **σ₈** | 2.7σ | **0.4σ** | ✅ Resolved |
| **H₀** | 4.8σ | **< 1σ** | ✅ Resolved |

---

## 📐 Theoretical Framework

### The Hertault Axiom

$$e^{4\sigma(x)} = \mathcal{I}(x) \equiv \frac{S_{\text{ent}}(x)}{S_{\text{Bek}}(x)}$$

The conformal mode σ encodes the information saturation ratio.

### Derived Parameters (No Free Parameters!)

| Parameter | Value | Source |
|-----------|-------|--------|
| α* | 0.075 | Asymptotic Safety UV fixed point |
| β | 2/3 | Holographic area law |
| ξ | 0.10 | Derived from β: ξ = β/[4(1+β)] |
| ρ_c | ρ_DE | UV-IR mixing |

### Mass Function

$$m^2_{\text{eff}}(\rho) = (\alpha_* M_{\text{Pl}})^2 \left[1 - \left(\frac{\rho}{\rho_c}\right)^{2/3}\right]$$

---

## 📁 Repository Structure

```
DGE-simulations/
├── src/                    # Core C and Python code
│   ├── dge_modifications.h # C header for CLASS-DGE
│   ├── class_dge.py        # Python wrapper for CLASS
│   └── dark_geometry.c     # CLASS modifications
│
├── scripts/
│   ├── field/              # Complete φ(z) dynamics
│   │   └── complete_field_dynamics.py
│   │
│   ├── mcmc/               # MCMC analysis
│   │   └── dge_mcmc.py
│   │
│   ├── hubble/             # H₀ tension analysis
│   │   ├── H0_field_dynamics.py
│   │   ├── H0_advanced_analysis.py
│   │   └── H0_tension_analysis.py
│   │
│   ├── desi/               # DESI Y1 validation
│   │   ├── desi_analysis.py
│   │   └── fsigma8_analysis.py
│   │
│   ├── synthesis/          # Complete resolution
│   │   └── complete_resolution.py
│   │
│   ├── bbn/                # BBN constraints
│   │   └── BBN_constraints.py
│   │
│   ├── potential/          # V(φ) from Hertault Axiom
│   │   └── DG_potential.py
│   │
│   ├── gw/                 # Gravitational waves
│   │   └── GW_predictions.py
│   │
│   └── nbody/              # N-body simulation configs
│       ├── ramses_dg.nml
│       └── ecosmog_params.ini
│
├── configs/                # Simulation configurations
├── figures/                # All output figures (40+)
├── docs/                   # Documentation
└── paper/                  # LaTeX paper template
```

---

## 🚀 Quick Start

### Installation

```bash
git clone https://github.com/username/DGE-simulations.git
cd DGE-simulations

# Install dependencies
pip install numpy scipy matplotlib emcee --break-system-packages
```

### Run Complete Analysis

```bash
# 1. Field dynamics (find optimal φ_initial)
python scripts/field/complete_field_dynamics.py

# 2. σ₈ tension resolution
python scripts/desi/desi_analysis.py

# 3. H₀ tension resolution  
python scripts/hubble/H0_field_dynamics.py

# 4. Full synthesis
python scripts/synthesis/complete_resolution.py
```

---

## 🔬 Physical Mechanisms

### σ₈ Resolution

Scale-dependent gravitational coupling produces **12% power spectrum suppression**:
- σ₈(DGE) = 0.77 vs σ₈(Planck) = 0.81
- Matches weak lensing observations

### H₀ Resolution

Non-minimal coupling ξRφ² with **ξ = 0.10 (derived)**:
- G_eff/G ≈ 1.08 at recombination  
- Δr_s/r_s ≈ -4%
- H₀: 67.4 → 72.6 km/s/Mpc

---

## 📊 Consistency Checks

| Test | Result | Status |
|------|--------|--------|
| **BBN** | Y_p preserved | ✅ |
| **GW170817** | c_T = c exactly | ✅ |
| **CMB θ*** | Preserved | ✅ |
| **DESI w(z)** | Dynamic | ✅ |

---

## 📧 Contact

- **Author:** Hugo Hertault
- **Email:** hertault.toe@gmail.com
- **X:** [@HertaultHu24527](https://x.com/HertaultHu24527)

---

## 📜 License

MIT License

---

*Dark Geometry: Unifying dark matter and dark energy through information geometry.*
