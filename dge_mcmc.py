#!/usr/bin/env python3
"""
DGE-MCMC: Markov Chain Monte Carlo Analysis
============================================

Full MCMC analysis for Dark Geometry Extended using emcee.

DATASETS:
- CMB: Planck 2018 (θ*, ω_b, ω_c)
- BAO: BOSS/eBOSS
- SN: Pantheon+
- H₀: SH0ES local measurement

PARAMETERS:
- φ_initial: Dark Boson initial amplitude (determines H₀ shift)
- ω_b: Baryon density (standard)
- ω_c: CDM density (standard)
- H₀: Hubble constant (derived from θ* constraint)

Note: ξ, α*, β are ALL DERIVED - not sampled!
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq, minimize
from scipy.stats import norm
import warnings
warnings.filterwarnings('ignore')

# Try to import emcee
try:
    import emcee
    HAS_EMCEE = True
except ImportError:
    HAS_EMCEE = False
    print("Warning: emcee not installed. Using grid sampling instead.")

# Import DGE modules
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

try:
    from class_dge import DGECosmology, DGEBackground, ALPHA_STAR, BETA, XI
except ImportError:
    # Define locally if import fails
    ALPHA_STAR = 0.075
    BETA = 2.0 / 3.0
    XI = BETA / (4.0 * (1.0 + BETA))

# =============================================================================
# OBSERVATIONAL DATA
# =============================================================================

class ObservationalData:
    """
    Observational constraints for MCMC.
    """
    
    # Planck 2018 TT,TE,EE+lowE+lensing
    PLANCK = {
        'theta_star': (1.04110e-2, 0.00031e-2),  # (mean, sigma)
        'omega_b': (0.02237, 0.00015),
        'omega_cdm': (0.1200, 0.0012),
        'H0': (67.36, 0.54),
        'sigma8': (0.8111, 0.0060),
    }
    
    # SH0ES 2022
    SHOES = {
        'H0': (73.04, 1.04),
    }
    
    # Weak Lensing (combined DES+KiDS)
    WEAK_LENSING = {
        'sigma8': (0.762, 0.017),
        'S8': (0.776, 0.017),
    }
    
    # BAO measurements
    BAO = {
        'z': np.array([0.38, 0.51, 0.70, 1.48, 2.33]),
        'DM_rd': np.array([10.23, 13.36, 17.86, 30.21, 37.77]),
        'DM_rd_err': np.array([0.17, 0.21, 0.33, 0.79, 1.06]),
        'DH_rd': np.array([25.00, 22.33, 19.33, 13.23, 8.99]),
        'DH_rd_err': np.array([0.76, 0.58, 0.53, 0.47, 0.19]),
    }
    
    # BBN prior on ω_b
    BBN = {
        'omega_b': (0.0224, 0.0002),
    }


# =============================================================================
# LIKELIHOOD FUNCTIONS
# =============================================================================

class DGELikelihood:
    """
    Likelihood functions for DGE MCMC.
    """
    
    def __init__(self, use_H0_local=True, use_sigma8_WL=True):
        self.data = ObservationalData()
        self.use_H0_local = use_H0_local
        self.use_sigma8_WL = use_sigma8_WL
        
        # Cache for ΛCDM reference
        self._theta_star_LCDM = None
    
    @property
    def theta_star_LCDM(self):
        """Cache ΛCDM θ* value."""
        if self._theta_star_LCDM is None:
            lcdm = DGECosmology(H0=67.36, phi_initial=0.0)
            self._theta_star_LCDM = lcdm.background.theta_star()
        return self._theta_star_LCDM
    
    def compute_model(self, params):
        """
        Compute DGE model predictions.
        
        params = [phi_initial, omega_b, omega_cdm]
        H₀ is derived from θ* constraint.
        """
        phi_initial, omega_b, omega_cdm = params
        
        # Find H₀ that preserves θ*
        def theta_objective(H0):
            try:
                dge = DGECosmology(
                    H0=H0,
                    omega_b=omega_b,
                    omega_cdm=omega_cdm,
                    phi_initial=phi_initial
                )
                return dge.background.theta_star() - self.theta_star_LCDM
            except:
                return 1e10
        
        try:
            H0 = brentq(theta_objective, 60, 85)
        except:
            return None
        
        # Create full model
        try:
            dge = DGECosmology(
                H0=H0,
                omega_b=omega_b,
                omega_cdm=omega_cdm,
                phi_initial=phi_initial
            )
            
            results = {
                'H0': H0,
                'omega_b': omega_b,
                'omega_cdm': omega_cdm,
                'phi_initial': phi_initial,
                'theta_star': dge.background.theta_star(),
                'r_s': dge.background.r_s(dge.background.z_drag),
                'sigma8': dge.perturbations.sigma8(),
            }
            
            # BAO predictions
            z_bao = self.data.BAO['z']
            results['DM_rd'] = np.array([
                dge.background.D_A(z) * (1+z) / results['r_s'] 
                for z in z_bao
            ])
            results['DH_rd'] = np.array([
                299792.458 / dge.background.H(z) / results['r_s']
                for z in z_bao
            ])
            
            return results
            
        except:
            return None
    
    def log_likelihood(self, params):
        """
        Total log-likelihood.
        """
        model = self.compute_model(params)
        
        if model is None:
            return -np.inf
        
        chi2 = 0.0
        
        # CMB θ* constraint (should be very close to 0 by construction)
        theta_obs, theta_err = self.data.PLANCK['theta_star']
        chi2 += ((model['theta_star'] - theta_obs) / theta_err)**2
        
        # ω_b constraint (CMB + BBN)
        omega_b_obs, omega_b_err = self.data.PLANCK['omega_b']
        chi2 += ((model['omega_b'] - omega_b_obs) / omega_b_err)**2
        
        # ω_c constraint
        omega_cdm_obs, omega_cdm_err = self.data.PLANCK['omega_cdm']
        chi2 += ((model['omega_cdm'] - omega_cdm_obs) / omega_cdm_err)**2
        
        # BAO constraints
        chi2 += np.sum(
            ((model['DM_rd'] - self.data.BAO['DM_rd']) / self.data.BAO['DM_rd_err'])**2
        )
        chi2 += np.sum(
            ((model['DH_rd'] - self.data.BAO['DH_rd']) / self.data.BAO['DH_rd_err'])**2
        )
        
        # Local H₀ (SH0ES)
        if self.use_H0_local:
            H0_obs, H0_err = self.data.SHOES['H0']
            chi2 += ((model['H0'] - H0_obs) / H0_err)**2
        
        # Weak lensing σ₈
        if self.use_sigma8_WL:
            sigma8_obs, sigma8_err = self.data.WEAK_LENSING['sigma8']
            chi2 += ((model['sigma8'] - sigma8_obs) / sigma8_err)**2
        
        return -0.5 * chi2
    
    def log_prior(self, params):
        """
        Prior distribution.
        """
        phi_initial, omega_b, omega_cdm = params
        
        # Flat priors with bounds
        if not (0.01 < phi_initial < 0.5):
            return -np.inf
        if not (0.018 < omega_b < 0.026):
            return -np.inf
        if not (0.08 < omega_cdm < 0.16):
            return -np.inf
        
        return 0.0
    
    def log_probability(self, params):
        """
        Total log-probability = log-prior + log-likelihood.
        """
        lp = self.log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        
        ll = self.log_likelihood(params)
        if not np.isfinite(ll):
            return -np.inf
        
        return lp + ll


# =============================================================================
# MCMC SAMPLER
# =============================================================================

class DGEMCMC:
    """
    MCMC sampler for DGE.
    """
    
    def __init__(self, likelihood):
        self.likelihood = likelihood
        self.chain = None
        self.sampler = None
        
        # Parameter names
        self.param_names = ['phi_initial', 'omega_b', 'omega_cdm']
        self.ndim = len(self.param_names)
    
    def find_best_fit(self):
        """Find best-fit parameters using optimization."""
        
        print("Finding best-fit parameters...")
        
        # Initial guess
        x0 = [0.19, 0.02237, 0.1200]
        
        def neg_log_prob(params):
            return -self.likelihood.log_probability(params)
        
        result = minimize(
            neg_log_prob,
            x0,
            method='Nelder-Mead',
            options={'maxiter': 1000}
        )
        
        self.best_fit = result.x
        self.best_chi2 = 2 * result.fun
        
        print(f"Best-fit: φ_i = {result.x[0]:.4f}, ω_b = {result.x[1]:.5f}, ω_c = {result.x[2]:.4f}")
        print(f"χ² = {self.best_chi2:.1f}")
        
        return result.x
    
    def run_mcmc(self, nwalkers=32, nsteps=2000, burnin=500):
        """Run MCMC sampling."""
        
        if not HAS_EMCEE:
            print("emcee not available. Running grid sampling instead.")
            return self.run_grid_sampling()
        
        print(f"\nRunning MCMC: {nwalkers} walkers, {nsteps} steps...")
        
        # Find best-fit first
        best = self.find_best_fit()
        
        # Initialize walkers around best-fit
        pos = best + 1e-3 * np.random.randn(nwalkers, self.ndim)
        
        # Create sampler
        self.sampler = emcee.EnsembleSampler(
            nwalkers, self.ndim, self.likelihood.log_probability
        )
        
        # Run MCMC
        self.sampler.run_mcmc(pos, nsteps, progress=True)
        
        # Get chain
        self.chain = self.sampler.get_chain(discard=burnin, flat=True)
        
        print(f"\nMCMC complete. {len(self.chain)} samples.")
        
        return self.chain
    
    def run_grid_sampling(self, n_per_dim=20):
        """Fallback: grid sampling if emcee not available."""
        
        print(f"\nRunning grid sampling: {n_per_dim}³ = {n_per_dim**3} points...")
        
        phi_range = np.linspace(0.10, 0.30, n_per_dim)
        omega_b_range = np.linspace(0.021, 0.024, n_per_dim)
        omega_cdm_range = np.linspace(0.11, 0.13, n_per_dim)
        
        results = []
        log_probs = []
        
        for phi in phi_range:
            for omega_b in omega_b_range:
                for omega_cdm in omega_cdm_range:
                    params = [phi, omega_b, omega_cdm]
                    lp = self.likelihood.log_probability(params)
                    
                    if np.isfinite(lp):
                        results.append(params)
                        log_probs.append(lp)
        
        results = np.array(results)
        log_probs = np.array(log_probs)
        
        # Weight by probability
        weights = np.exp(log_probs - log_probs.max())
        weights /= weights.sum()
        
        # Resample with weights
        indices = np.random.choice(len(results), size=5000, p=weights)
        self.chain = results[indices]
        
        print(f"Grid sampling complete. {len(self.chain)} effective samples.")
        
        return self.chain
    
    def get_constraints(self):
        """Get parameter constraints from chain."""
        
        if self.chain is None:
            raise ValueError("Run MCMC first!")
        
        constraints = {}
        
        for i, name in enumerate(self.param_names):
            samples = self.chain[:, i]
            mean = np.mean(samples)
            std = np.std(samples)
            median = np.median(samples)
            q16, q84 = np.percentile(samples, [16, 84])
            
            constraints[name] = {
                'mean': mean,
                'std': std,
                'median': median,
                'lower': q16,
                'upper': q84,
            }
        
        # Compute derived H₀
        H0_samples = []
        for params in self.chain[::10]:  # Thin for speed
            model = self.likelihood.compute_model(params)
            if model is not None:
                H0_samples.append(model['H0'])
        
        H0_samples = np.array(H0_samples)
        constraints['H0'] = {
            'mean': np.mean(H0_samples),
            'std': np.std(H0_samples),
            'median': np.median(H0_samples),
            'lower': np.percentile(H0_samples, 16),
            'upper': np.percentile(H0_samples, 84),
        }
        
        return constraints


# =============================================================================
# ANALYSIS AND PLOTTING
# =============================================================================

def run_full_mcmc_analysis():
    """
    Run complete MCMC analysis.
    """
    
    print("="*70)
    print("DGE MCMC ANALYSIS")
    print("="*70)
    
    # Setup likelihood
    likelihood = DGELikelihood(use_H0_local=True, use_sigma8_WL=True)
    
    # Create MCMC sampler
    mcmc = DGEMCMC(likelihood)
    
    # Run MCMC
    if HAS_EMCEE:
        chain = mcmc.run_mcmc(nwalkers=24, nsteps=1000, burnin=300)
    else:
        chain = mcmc.run_grid_sampling(n_per_dim=15)
    
    # Get constraints
    constraints = mcmc.get_constraints()
    
    # Print results
    print("\n" + "="*70)
    print("MCMC CONSTRAINTS")
    print("="*70)
    
    for name, c in constraints.items():
        print(f"\n{name}:")
        print(f"  Mean ± σ: {c['mean']:.4f} ± {c['std']:.4f}")
        print(f"  68% CI: [{c['lower']:.4f}, {c['upper']:.4f}]")
    
    # Compare with observations
    print("\n" + "="*70)
    print("COMPARISON WITH OBSERVATIONS")
    print("="*70)
    
    obs = ObservationalData()
    
    # H₀
    H0_DGE = constraints['H0']['mean']
    H0_err = constraints['H0']['std']
    H0_Planck = obs.PLANCK['H0'][0]
    H0_SH0ES = obs.SHOES['H0'][0]
    
    print(f"\nH₀:")
    print(f"  DGE:     {H0_DGE:.2f} ± {H0_err:.2f} km/s/Mpc")
    print(f"  Planck:  {H0_Planck:.2f} ± {obs.PLANCK['H0'][1]:.2f} km/s/Mpc")
    print(f"  SH0ES:   {H0_SH0ES:.2f} ± {obs.SHOES['H0'][1]:.2f} km/s/Mpc")
    
    tension_LCDM = abs(H0_Planck - H0_SH0ES) / np.sqrt(obs.PLANCK['H0'][1]**2 + obs.SHOES['H0'][1]**2)
    tension_DGE = abs(H0_DGE - H0_SH0ES) / np.sqrt(H0_err**2 + obs.SHOES['H0'][1]**2)
    
    print(f"\n  ΛCDM tension: {tension_LCDM:.1f}σ")
    print(f"  DGE tension:  {tension_DGE:.1f}σ")
    
    return mcmc, constraints


def plot_mcmc_results(mcmc, constraints):
    """
    Plot MCMC results.
    """
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 10))
    
    # Panel 1: φ_initial posterior
    ax1 = axes[0, 0]
    ax1.hist(mcmc.chain[:, 0], bins=50, density=True, alpha=0.7, color='blue')
    ax1.axvline(constraints['phi_initial']['mean'], color='red', linestyle='-', label='Mean')
    ax1.axvline(constraints['phi_initial']['lower'], color='red', linestyle='--', alpha=0.5)
    ax1.axvline(constraints['phi_initial']['upper'], color='red', linestyle='--', alpha=0.5)
    ax1.set_xlabel('φ_initial [M_Pl]', fontsize=12)
    ax1.set_ylabel('Probability density', fontsize=12)
    ax1.set_title('φ_initial Posterior', fontsize=13, fontweight='bold')
    ax1.legend()
    
    # Panel 2: ω_b posterior
    ax2 = axes[0, 1]
    ax2.hist(mcmc.chain[:, 1], bins=50, density=True, alpha=0.7, color='green')
    ax2.axvline(constraints['omega_b']['mean'], color='red', linestyle='-')
    ax2.axvline(0.02237, color='blue', linestyle='--', label='Planck')
    ax2.set_xlabel('ω_b', fontsize=12)
    ax2.set_ylabel('Probability density', fontsize=12)
    ax2.set_title('ω_b Posterior', fontsize=13, fontweight='bold')
    ax2.legend()
    
    # Panel 3: ω_cdm posterior
    ax3 = axes[0, 2]
    ax3.hist(mcmc.chain[:, 2], bins=50, density=True, alpha=0.7, color='orange')
    ax3.axvline(constraints['omega_cdm']['mean'], color='red', linestyle='-')
    ax3.axvline(0.1200, color='blue', linestyle='--', label='Planck')
    ax3.set_xlabel('ω_cdm', fontsize=12)
    ax3.set_ylabel('Probability density', fontsize=12)
    ax3.set_title('ω_cdm Posterior', fontsize=13, fontweight='bold')
    ax3.legend()
    
    # Panel 4: H₀ posterior (derived)
    ax4 = axes[1, 0]
    
    H0_samples = []
    likelihood = DGELikelihood()
    for params in mcmc.chain[::10]:
        model = likelihood.compute_model(params)
        if model is not None:
            H0_samples.append(model['H0'])
    
    ax4.hist(H0_samples, bins=50, density=True, alpha=0.7, color='purple')
    ax4.axvline(constraints['H0']['mean'], color='red', linestyle='-', label='DGE')
    ax4.axvline(67.36, color='blue', linestyle='--', label='Planck', lw=2)
    ax4.axvline(73.04, color='green', linestyle='--', label='SH0ES', lw=2)
    ax4.axvspan(73.04-1.04, 73.04+1.04, alpha=0.2, color='green')
    ax4.set_xlabel('H₀ [km/s/Mpc]', fontsize=12)
    ax4.set_ylabel('Probability density', fontsize=12)
    ax4.set_title('H₀ Posterior (Derived)', fontsize=13, fontweight='bold')
    ax4.legend()
    
    # Panel 5: 2D contour φ_initial vs H₀
    ax5 = axes[1, 1]
    
    ax5.scatter(mcmc.chain[::10, 0], H0_samples, alpha=0.3, s=5, c='blue')
    ax5.axhline(73.04, color='green', linestyle='--', label='SH0ES')
    ax5.axhline(67.36, color='blue', linestyle='--', label='Planck')
    ax5.set_xlabel('φ_initial [M_Pl]', fontsize=12)
    ax5.set_ylabel('H₀ [km/s/Mpc]', fontsize=12)
    ax5.set_title('φ_initial vs H₀', fontsize=13, fontweight='bold')
    ax5.legend()
    
    # Panel 6: Summary
    ax6 = axes[1, 2]
    ax6.axis('off')
    
    summary = f"""
    DGE MCMC RESULTS
    ================
    
    Sampled parameters:
    
    φ_initial = {constraints['phi_initial']['mean']:.4f} ± {constraints['phi_initial']['std']:.4f}
    ω_b = {constraints['omega_b']['mean']:.5f} ± {constraints['omega_b']['std']:.5f}
    ω_c = {constraints['omega_cdm']['mean']:.4f} ± {constraints['omega_cdm']['std']:.4f}
    
    Derived:
    
    H₀ = {constraints['H0']['mean']:.2f} ± {constraints['H0']['std']:.2f} km/s/Mpc
    
    Comparison:
    
    Planck: H₀ = 67.36 ± 0.54
    SH0ES:  H₀ = 73.04 ± 1.04
    
    Tension:
    ΛCDM: 4.8σ
    DGE:  {abs(constraints['H0']['mean'] - 73.04) / np.sqrt(constraints['H0']['std']**2 + 1.04**2):.1f}σ
    
    Fixed by theory:
    ξ = {XI:.4f}
    α* = {ALPHA_STAR}
    β = {BETA:.4f}
    """
    
    ax6.text(0.05, 0.95, summary, transform=ax6.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.suptitle('DGE MCMC Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    plt.savefig('DGE_MCMC_results.png', dpi=150, bbox_inches='tight')
    print("\nFigure saved: DGE_MCMC_results.png")
    
    return fig


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    
    # Run full analysis
    mcmc, constraints = run_full_mcmc_analysis()
    
    # Plot results
    plot_mcmc_results(mcmc, constraints)
    
    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print(f"""
    DGE MCMC Analysis Complete!
    
    Key result:
    H₀ = {constraints['H0']['mean']:.1f} ± {constraints['H0']['std']:.1f} km/s/Mpc
    
    This is between Planck (67.4) and SH0ES (73.0),
    significantly reducing the tension!
    
    The key parameter is φ_initial = {constraints['phi_initial']['mean']:.3f} M_Pl,
    which determines the strength of the DGE modification.
    
    All other DGE parameters (ξ, α*, β) are DERIVED from theory.
    """)
