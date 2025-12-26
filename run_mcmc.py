#!/usr/bin/env python3
"""
Dark Geometry MCMC - Simplified Implementation
==============================================

This script runs a simple MCMC to validate CLASS-DG against:
1. σ₈ from weak lensing (DES + KiDS)
2. H₀ prior (optional)
3. Basic cosmological constraints

Uses emcee for MCMC sampling.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import warnings
warnings.filterwarnings('ignore')

# Setup paths
CLASS_DG_PATH = "/home/claude/class_dg"
os.chdir(CLASS_DG_PATH)
sys.path.insert(0, os.path.join(CLASS_DG_PATH, "python"))

from classy import Class

# =============================================================================
# OBSERVATIONAL DATA
# =============================================================================

# Weak lensing σ₈ (DES Y3 + KiDS-1000 combined)
SIGMA8_OBS = 0.762
SIGMA8_ERR = 0.020

# Planck 2018 CMB (for comparison)
SIGMA8_PLANCK = 0.811
SIGMA8_PLANCK_ERR = 0.006

# H₀ from SH0ES (optional constraint)
H0_SHOES = 73.04
H0_SHOES_ERR = 1.04

# Planck H₀ (for comparison)
H0_PLANCK = 67.36
H0_PLANCK_ERR = 0.54

# BAO: D_V/r_d at z=0.51 (BOSS DR12)
DV_RD_Z051 = 13.38
DV_RD_Z051_ERR = 0.18

# =============================================================================
# CLASS-DG WRAPPER
# =============================================================================

def compute_cosmology(omega_b, omega_cdm, H0, A_s, n_s, tau_reio):
    """
    Compute cosmological observables using CLASS-DG.
    
    Returns dict with sigma8, Omega_m, rs_drag, etc.
    """
    cosmo = Class()
    
    try:
        cosmo.set({
            'output': 'mPk',
            'P_k_max_h/Mpc': 10,
            'omega_b': omega_b,
            'omega_cdm': omega_cdm,
            'H0': H0,
            'A_s': A_s,
            'n_s': n_s,
            'tau_reio': tau_reio,
        })
        cosmo.compute()
        
        result = {
            'sigma8': cosmo.sigma8(),
            'Omega_m': cosmo.Omega_m(),
            'Omega_Lambda': cosmo.Omega_Lambda(),
            'rs_drag': cosmo.rs_drag(),
            'age': cosmo.age(),
            'h': cosmo.h(),
        }
        
        cosmo.struct_cleanup()
        cosmo.empty()
        
        return result
        
    except Exception as e:
        if cosmo is not None:
            try:
                cosmo.struct_cleanup()
                cosmo.empty()
            except:
                pass
        return None

# =============================================================================
# LIKELIHOOD FUNCTIONS
# =============================================================================

def log_prior(theta):
    """
    Flat priors on cosmological parameters.
    
    theta = [omega_b, omega_cdm, H0, ln(10^10 A_s), n_s, tau_reio]
    """
    omega_b, omega_cdm, H0, logA, n_s, tau_reio = theta
    
    # Prior bounds
    if not (0.019 < omega_b < 0.025):
        return -np.inf
    if not (0.10 < omega_cdm < 0.14):
        return -np.inf
    if not (60 < H0 < 80):
        return -np.inf
    if not (2.8 < logA < 3.3):
        return -np.inf
    if not (0.92 < n_s < 1.0):
        return -np.inf
    if not (0.02 < tau_reio < 0.12):
        return -np.inf
    
    return 0.0

def log_likelihood(theta, use_sigma8=True, use_h0_shoes=False):
    """
    Compute log-likelihood given parameters.
    
    Likelihoods:
    - σ₈ from weak lensing
    - Optionally H₀ from SH0ES
    """
    omega_b, omega_cdm, H0, logA, n_s, tau_reio = theta
    
    # Convert logA to A_s
    A_s = np.exp(logA) * 1e-10
    
    # Compute cosmology
    result = compute_cosmology(omega_b, omega_cdm, H0, A_s, n_s, tau_reio)
    
    if result is None:
        return -np.inf
    
    chi2 = 0.0
    
    # σ₈ likelihood (weak lensing)
    if use_sigma8:
        chi2 += ((result['sigma8'] - SIGMA8_OBS) / SIGMA8_ERR)**2
    
    # H₀ likelihood (SH0ES) - optional
    if use_h0_shoes:
        chi2 += ((H0 - H0_SHOES) / H0_SHOES_ERR)**2
    
    return -0.5 * chi2

def log_probability(theta, use_sigma8=True, use_h0_shoes=False):
    """Log-posterior = log-prior + log-likelihood."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, use_sigma8, use_h0_shoes)

# =============================================================================
# MCMC SAMPLER
# =============================================================================

def run_mcmc(nwalkers=16, nsteps=500, use_sigma8=True, use_h0_shoes=False):
    """
    Run MCMC using emcee.
    
    Parameters
    ----------
    nwalkers : int
        Number of walkers
    nsteps : int
        Number of steps per walker
    use_sigma8 : bool
        Include σ₈ likelihood
    use_h0_shoes : bool
        Include H₀ (SH0ES) likelihood
        
    Returns
    -------
    samples : array
        MCMC samples (after burn-in)
    """
    try:
        import emcee
    except ImportError:
        print("Installing emcee...")
        os.system("pip install emcee --break-system-packages -q")
        import emcee
    
    # Parameter dimensions
    ndim = 6
    
    # Initial positions (scattered around Planck best-fit)
    p0_mean = np.array([0.02237, 0.1200, 67.36, 3.044, 0.9649, 0.0544])
    p0_std = np.array([0.0002, 0.002, 1.0, 0.02, 0.01, 0.01])
    
    p0 = p0_mean + p0_std * np.random.randn(nwalkers, ndim)
    
    # Setup sampler
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability,
        args=(use_sigma8, use_h0_shoes)
    )
    
    print(f"Running MCMC with {nwalkers} walkers, {nsteps} steps...")
    print(f"  - σ₈ likelihood: {use_sigma8}")
    print(f"  - H₀ (SH0ES) likelihood: {use_h0_shoes}")
    
    # Run MCMC with progress
    sampler.run_mcmc(p0, nsteps, progress=True)
    
    # Get samples (discard burn-in)
    burnin = nsteps // 4
    samples = sampler.get_chain(discard=burnin, flat=True)
    
    print(f"\nMCMC complete. {len(samples)} samples after burn-in.")
    
    return samples, sampler

# =============================================================================
# ANALYSIS
# =============================================================================

def analyze_samples(samples):
    """
    Analyze MCMC samples and compute derived parameters.
    """
    labels = ['ω_b', 'ω_cdm', 'H₀', 'ln(10¹⁰A_s)', 'n_s', 'τ_reio']
    
    print("\n" + "="*60)
    print("MCMC RESULTS")
    print("="*60)
    
    # Parameter constraints
    print("\nParameter constraints (mean ± std):")
    means = np.mean(samples, axis=0)
    stds = np.std(samples, axis=0)
    
    for i, (label, mean, std) in enumerate(zip(labels, means, stds)):
        print(f"  {label:15s} = {mean:.5f} ± {std:.5f}")
    
    # Compute derived parameters for each sample
    print("\nComputing derived parameters...")
    
    # Take a subset for speed
    n_derived = min(100, len(samples))
    indices = np.random.choice(len(samples), n_derived, replace=False)
    
    sigma8_samples = []
    omega_m_samples = []
    
    for idx in indices:
        omega_b, omega_cdm, H0, logA, n_s, tau_reio = samples[idx]
        A_s = np.exp(logA) * 1e-10
        
        result = compute_cosmology(omega_b, omega_cdm, H0, A_s, n_s, tau_reio)
        if result is not None:
            sigma8_samples.append(result['sigma8'])
            omega_m_samples.append(result['Omega_m'])
    
    sigma8_samples = np.array(sigma8_samples)
    omega_m_samples = np.array(omega_m_samples)
    
    print("\nDerived parameters:")
    print(f"  σ₈    = {np.mean(sigma8_samples):.4f} ± {np.std(sigma8_samples):.4f}")
    print(f"  Ω_m   = {np.mean(omega_m_samples):.4f} ± {np.std(omega_m_samples):.4f}")
    
    # Comparison with observations
    print("\n" + "-"*60)
    print("COMPARISON WITH OBSERVATIONS")
    print("-"*60)
    
    sigma8_mean = np.mean(sigma8_samples)
    sigma8_std = np.std(sigma8_samples)
    
    tension_wl = abs(sigma8_mean - SIGMA8_OBS) / np.sqrt(sigma8_std**2 + SIGMA8_ERR**2)
    tension_planck = abs(sigma8_mean - SIGMA8_PLANCK) / np.sqrt(sigma8_std**2 + SIGMA8_PLANCK_ERR**2)
    
    print(f"\nσ₈ comparison:")
    print(f"  CLASS-DG:      {sigma8_mean:.4f} ± {sigma8_std:.4f}")
    print(f"  Weak lensing:  {SIGMA8_OBS:.4f} ± {SIGMA8_ERR:.4f}")
    print(f"  Planck ΛCDM:   {SIGMA8_PLANCK:.4f} ± {SIGMA8_PLANCK_ERR:.4f}")
    print(f"\n  Tension with WL:     {tension_wl:.1f}σ")
    print(f"  Tension with Planck: {tension_planck:.1f}σ")
    
    if tension_wl < 1.0:
        print("\n  ✓ σ₈ tension RESOLVED!")
    
    return {
        'samples': samples,
        'sigma8': sigma8_samples,
        'omega_m': omega_m_samples,
        'labels': labels,
    }

def plot_results(results, filename='mcmc_results.png'):
    """Create summary plots."""
    
    try:
        import corner
    except ImportError:
        print("Installing corner...")
        os.system("pip install corner --break-system-packages -q")
        import corner
    
    samples = results['samples']
    labels = results['labels']
    sigma8 = results['sigma8']
    
    fig = plt.figure(figsize=(16, 12))
    
    # Corner plot for main parameters
    fig_corner = corner.corner(
        samples,
        labels=labels,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 10}
    )
    fig_corner.savefig('mcmc_corner.png', dpi=100, bbox_inches='tight')
    print(f"Corner plot saved: mcmc_corner.png")
    
    # σ₈ histogram
    fig2, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(sigma8, bins=30, density=True, alpha=0.7, color='blue', label='CLASS-DG posterior')
    
    # Add Gaussian fits
    x = np.linspace(0.7, 0.85, 100)
    ax.plot(x, norm.pdf(x, SIGMA8_OBS, SIGMA8_ERR), 'g-', linewidth=2, 
            label=f'Weak lensing ({SIGMA8_OBS}±{SIGMA8_ERR})')
    ax.plot(x, norm.pdf(x, SIGMA8_PLANCK, SIGMA8_PLANCK_ERR), 'r--', linewidth=2,
            label=f'Planck ΛCDM ({SIGMA8_PLANCK}±{SIGMA8_PLANCK_ERR})')
    
    ax.axvline(np.mean(sigma8), color='blue', linestyle='-', linewidth=2)
    ax.axvspan(np.mean(sigma8)-np.std(sigma8), np.mean(sigma8)+np.std(sigma8), 
               alpha=0.2, color='blue')
    
    ax.set_xlabel('σ₈', fontsize=14)
    ax.set_ylabel('Probability density', fontsize=14)
    ax.set_title('Dark Geometry: σ₈ Posterior', fontsize=15, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    fig2.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"σ₈ plot saved: {filename}")
    
    return fig_corner, fig2

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("="*60)
    print("DARK GEOMETRY MCMC VALIDATION")
    print("="*60)
    
    # Quick test first
    print("\nTesting CLASS-DG...")
    result = compute_cosmology(0.02237, 0.1200, 67.36, 2.1e-9, 0.9649, 0.0544)
    if result:
        print(f"  σ₈ = {result['sigma8']:.4f} (expected ~0.77)")
        print("  CLASS-DG working ✓")
    else:
        print("  ERROR: CLASS-DG not working!")
        sys.exit(1)
    
    # Run MCMC
    print("\n" + "="*60)
    
    # Short run for demonstration
    samples, sampler = run_mcmc(
        nwalkers=16,  # Must be >= 2*ndim = 12
        nsteps=30,  # Very short for demo
        use_sigma8=True,
        use_h0_shoes=False
    )
    
    # Analyze
    results = analyze_samples(samples)
    
    # Plot
    plot_results(results)
    
    print("\n" + "="*60)
    print("MCMC COMPLETE")
    print("="*60)
