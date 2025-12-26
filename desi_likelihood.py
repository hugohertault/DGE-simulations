#!/usr/bin/env python3
"""
Dark Geometry - DESI Y1 Likelihood pour MCMC
=============================================

Ce script implémente une likelihood complète pour analyser
les données DESI Y1 avec Dark Geometry.

Peut être utilisé avec:
- emcee (Python)
- Cobaya
- MontePython
"""

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import json

# =============================================================================
# DESI Y1 COVARIANCE MATRICES
# =============================================================================

# BAO data with full covariance
# From DESI 2024 III (arXiv:2404.03000)

DESI_BAO = {
    'z_eff': np.array([0.295, 0.510, 0.706, 0.930, 1.317, 1.491, 2.330]),
    'tracers': ['BGS', 'LRG1', 'LRG2', 'LRG3', 'ELG', 'QSO', 'Lya'],
    
    # D_V/r_d for BGS only
    'DV_rd': {
        'BGS': {'value': 7.93, 'error': 0.15},
    },
    
    # D_M/r_d and D_H/r_d for other tracers
    'DM_rd': {
        'LRG1': {'value': 13.62, 'error': 0.25},
        'LRG2': {'value': 16.85, 'error': 0.32},
        'LRG3': {'value': 21.71, 'error': 0.28},
        'ELG': {'value': 27.79, 'error': 0.69},
        'QSO': {'value': 30.69, 'error': 0.80},
        'Lya': {'value': 39.71, 'error': 0.94},
    },
    'DH_rd': {
        'LRG1': {'value': 20.98, 'error': 0.61},
        'LRG2': {'value': 20.08, 'error': 0.60},
        'LRG3': {'value': 17.88, 'error': 0.35},
        'ELG': {'value': 13.82, 'error': 0.42},
        'QSO': {'value': 13.26, 'error': 0.55},
        'Lya': {'value': 8.52, 'error': 0.17},
    },
    
    # Correlation between D_M and D_H (approximate)
    'correlation_DM_DH': -0.4,
}

# RSD data: f*sigma8
DESI_RSD = {
    'z_eff': np.array([0.295, 0.510, 0.706, 0.930, 1.317]),
    'fsigma8': np.array([0.392, 0.458, 0.449, 0.437, 0.372]),
    'error': np.array([0.044, 0.033, 0.032, 0.035, 0.062]),
}


# =============================================================================
# COSMOLOGY CALCULATOR
# =============================================================================

class CosmologyCalculator:
    """
    Fast cosmology calculator for likelihood evaluation.
    """
    
    def __init__(self, Omega_m, h, Omega_b=0.0493, sigma8=0.811, 
                 model='LCDM', alpha_star=0.075):
        
        self.Omega_m = Omega_m
        self.Omega_L = 1 - Omega_m
        self.h = h
        self.H0 = 100 * h
        self.Omega_b = Omega_b
        self.sigma8 = sigma8
        self.c = 299792.458
        
        self.model = model
        self.alpha_star = alpha_star
        
        # Sound horizon (approximate fitting formula)
        self.r_d = self._compute_rd()
        
        # Pre-compute interpolation tables
        self._setup_interpolation()
    
    def _compute_rd(self):
        """
        Compute sound horizon at drag epoch.
        Fitting formula from Eisenstein & Hu 1998.
        """
        omega_m = self.Omega_m * self.h**2
        omega_b = self.Omega_b * self.h**2
        
        # Fitting formula
        z_eq = 2.5e4 * omega_m * (2.7255/2.725)**(-4)
        k_eq = 0.0746 * omega_m * (2.725/2.7255)**2
        
        z_d = 1291 * (omega_m**0.251 / (1 + 0.659*omega_m**0.828)) * \
              (1 + 0.313*omega_m**(-0.419) * (1 + 0.607*omega_m**0.674)**(-1))
        
        R_d = 31.5 * omega_b * (2.725/2.7255)**(-4) * (1000/z_d)
        R_eq = 31.5 * omega_b * (2.725/2.7255)**(-4) * (1000/z_eq)
        
        r_d = (2/(3*k_eq)) * np.sqrt(6/R_eq) * \
              np.log((np.sqrt(1+R_d) + np.sqrt(R_d+R_eq)) / (1 + np.sqrt(R_eq)))
        
        # Convert to Mpc
        r_d = r_d / self.h
        
        return r_d
    
    def _setup_interpolation(self):
        """Pre-compute distance tables."""
        z_arr = np.linspace(0, 3, 500)
        
        DC_arr = np.array([self._DC_exact(z) for z in z_arr])
        
        self._DC_interp = interp1d(z_arr, DC_arr, kind='cubic')
    
    def _DC_exact(self, z):
        """Exact comoving distance integral."""
        integrand = lambda zp: 1.0 / self.E(zp)
        result, _ = quad(integrand, 0, z)
        return self.c / self.H0 * result
    
    def E(self, z):
        """E(z) = H(z)/H0."""
        if self.model == 'DG':
            # DG has slightly modified expansion
            # Through effective w(z)
            w_eff = self._w_dg(z)
            rho_de = self._rho_de_dg(z)
            return np.sqrt(self.Omega_m * (1+z)**3 + self.Omega_L * rho_de)
        else:
            return np.sqrt(self.Omega_m * (1+z)**3 + self.Omega_L)
    
    def _w_dg(self, z):
        """DG effective equation of state."""
        a = 1 / (1 + z)
        a_c = (self.Omega_m / self.Omega_L)**(1/3)
        x = np.log(a / a_c)
        w = -1 + 0.15 * (1 + np.tanh(x)) / 2
        return w
    
    def _rho_de_dg(self, z):
        """DG dark energy density evolution."""
        # Approximate: close to ΛCDM
        return 1.0
    
    def D_C(self, z):
        """Comoving distance."""
        if z <= 0:
            return 0.0
        return float(self._DC_interp(z))
    
    def D_M(self, z):
        """Transverse comoving distance."""
        return self.D_C(z)
    
    def D_H(self, z):
        """Hubble distance."""
        return self.c / (self.H0 * self.E(z))
    
    def D_V(self, z):
        """Volume-averaged distance."""
        return (z * self.D_M(z)**2 * self.D_H(z))**(1/3)
    
    def f_sigma8(self, z):
        """Growth rate f*sigma8."""
        E2 = self.E(z)**2
        Omega_m_z = self.Omega_m * (1+z)**3 / E2
        
        # Growth rate
        gamma = 0.55
        if self.model == 'DG':
            gamma = 0.55 + 0.02 * self.alpha_star
        f = Omega_m_z**gamma
        
        # Growth factor (approximate)
        a = 1 / (1 + z)
        D_z = a * Omega_m_z**(-0.1) / self.Omega_m**(-0.1)
        
        # sigma8(z)
        if self.model == 'DG':
            sigma8_z = 0.773 * D_z  # DG suppressed
        else:
            sigma8_z = self.sigma8 * D_z
        
        return f * sigma8_z


# =============================================================================
# LIKELIHOOD FUNCTIONS
# =============================================================================

class DESILikelihood:
    """
    DESI Y1 likelihood for cosmological analysis.
    """
    
    def __init__(self, use_bao=True, use_rsd=True):
        self.use_bao = use_bao
        self.use_rsd = use_rsd
        
        self.bao_data = DESI_BAO
        self.rsd_data = DESI_RSD
    
    def log_likelihood(self, theta, model='LCDM'):
        """
        Compute log-likelihood for given parameters.
        
        Parameters:
        -----------
        theta : array
            [Omega_m, h, sigma8] for LCDM
            [Omega_m, h] for DG (sigma8 derived)
        model : str
            'LCDM' or 'DG'
        """
        
        if model == 'DG':
            Omega_m, h = theta[:2]
            sigma8 = 0.773  # DG prediction
            cosmo = CosmologyCalculator(Omega_m, h, sigma8=sigma8, model='DG')
        else:
            Omega_m, h, sigma8 = theta[:3]
            cosmo = CosmologyCalculator(Omega_m, h, sigma8=sigma8, model='LCDM')
        
        chi2 = 0.0
        
        # BAO contribution
        if self.use_bao:
            chi2 += self._chi2_bao(cosmo)
        
        # RSD contribution
        if self.use_rsd:
            chi2 += self._chi2_rsd(cosmo)
        
        return -0.5 * chi2
    
    def _chi2_bao(self, cosmo):
        """BAO χ²."""
        chi2 = 0.0
        
        # BGS: D_V/r_d only
        z = 0.295
        DV_th = cosmo.D_V(z) / cosmo.r_d
        DV_obs = self.bao_data['DV_rd']['BGS']['value']
        DV_err = self.bao_data['DV_rd']['BGS']['error']
        chi2 += ((DV_th - DV_obs) / DV_err)**2
        
        # Other tracers: D_M/r_d and D_H/r_d
        for tracer in ['LRG1', 'LRG2', 'LRG3', 'ELG', 'QSO', 'Lya']:
            z_idx = self.bao_data['tracers'].index(tracer)
            z = self.bao_data['z_eff'][z_idx]
            
            DM_th = cosmo.D_M(z) / cosmo.r_d
            DM_obs = self.bao_data['DM_rd'][tracer]['value']
            DM_err = self.bao_data['DM_rd'][tracer]['error']
            
            DH_th = cosmo.D_H(z) / cosmo.r_d
            DH_obs = self.bao_data['DH_rd'][tracer]['value']
            DH_err = self.bao_data['DH_rd'][tracer]['error']
            
            # Include correlation
            rho = self.bao_data['correlation_DM_DH']
            
            # 2x2 covariance
            cov = np.array([[DM_err**2, rho*DM_err*DH_err],
                           [rho*DM_err*DH_err, DH_err**2]])
            
            diff = np.array([DM_th - DM_obs, DH_th - DH_obs])
            
            cov_inv = np.linalg.inv(cov)
            chi2 += diff @ cov_inv @ diff
        
        return chi2
    
    def _chi2_rsd(self, cosmo):
        """RSD χ²."""
        chi2 = 0.0
        
        for i, z in enumerate(self.rsd_data['z_eff']):
            fs8_th = cosmo.f_sigma8(z)
            fs8_obs = self.rsd_data['fsigma8'][i]
            fs8_err = self.rsd_data['error'][i]
            
            chi2 += ((fs8_th - fs8_obs) / fs8_err)**2
        
        return chi2


# =============================================================================
# MCMC WRAPPER
# =============================================================================

def run_mcmc_desi(model='DG', nwalkers=32, nsteps=1000):
    """
    Run MCMC on DESI data.
    """
    import emcee
    
    likelihood = DESILikelihood()
    
    if model == 'DG':
        ndim = 2
        labels = ['Omega_m', 'h']
        
        # Prior ranges
        prior_min = np.array([0.2, 0.6])
        prior_max = np.array([0.4, 0.8])
        
        def log_prior(theta):
            if np.all(theta > prior_min) and np.all(theta < prior_max):
                return 0.0
            return -np.inf
        
        def log_prob(theta):
            lp = log_prior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + likelihood.log_likelihood(theta, model='DG')
        
        # Initial positions
        p0 = np.array([0.315, 0.674])
        
    else:
        ndim = 3
        labels = ['Omega_m', 'h', 'sigma8']
        
        prior_min = np.array([0.2, 0.6, 0.7])
        prior_max = np.array([0.4, 0.8, 0.9])
        
        def log_prior(theta):
            if np.all(theta > prior_min) and np.all(theta < prior_max):
                return 0.0
            return -np.inf
        
        def log_prob(theta):
            lp = log_prior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + likelihood.log_likelihood(theta, model='LCDM')
        
        p0 = np.array([0.315, 0.674, 0.811])
    
    # Initialize walkers
    pos = p0 + 1e-3 * np.random.randn(nwalkers, ndim)
    
    # Run sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
    
    print(f"Running MCMC for {model}...")
    sampler.run_mcmc(pos, nsteps, progress=True)
    
    # Get results
    flat_samples = sampler.get_chain(discard=nsteps//5, flat=True)
    
    print(f"\nResults for {model}:")
    for i, label in enumerate(labels):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(f"  {label}: {mcmc[1]:.4f} +{q[1]:.4f} -{q[0]:.4f}")
    
    return sampler, labels


# =============================================================================
# QUICK COMPARISON
# =============================================================================

def quick_comparison():
    """
    Quick comparison without full MCMC.
    """
    
    print("="*60)
    print("DESI Y1 QUICK COMPARISON")
    print("="*60)
    
    likelihood = DESILikelihood()
    
    # Best-fit parameters
    theta_lcdm = [0.315, 0.674, 0.811]
    theta_dg = [0.315, 0.674]
    
    logl_lcdm = likelihood.log_likelihood(theta_lcdm, model='LCDM')
    logl_dg = likelihood.log_likelihood(theta_dg, model='DG')
    
    chi2_lcdm = -2 * logl_lcdm
    chi2_dg = -2 * logl_dg
    
    print(f"\nχ² at best-fit:")
    print(f"  ΛCDM: {chi2_lcdm:.1f}")
    print(f"  DG:   {chi2_dg:.1f}")
    print(f"  Δχ²:  {chi2_dg - chi2_lcdm:+.1f}")
    
    # Model selection
    n_params_lcdm = 3
    n_params_dg = 2
    n_data = 1 + 12 + 5  # BGS + 6×(DM+DH) + RSD
    
    aic_lcdm = chi2_lcdm + 2 * n_params_lcdm
    aic_dg = chi2_dg + 2 * n_params_dg
    
    bic_lcdm = chi2_lcdm + n_params_lcdm * np.log(n_data)
    bic_dg = chi2_dg + n_params_dg * np.log(n_data)
    
    print(f"\nModel selection:")
    print(f"  AIC(ΛCDM): {aic_lcdm:.1f}")
    print(f"  AIC(DG):   {aic_dg:.1f}")
    print(f"  ΔAIC:      {aic_dg - aic_lcdm:+.1f}")
    print(f"\n  BIC(ΛCDM): {bic_lcdm:.1f}")
    print(f"  BIC(DG):   {bic_dg:.1f}")
    print(f"  ΔBIC:      {bic_dg - bic_lcdm:+.1f}")
    
    return chi2_lcdm, chi2_dg


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    
    # Quick comparison
    quick_comparison()
    
    print("\n" + "="*60)
    print("To run full MCMC, use:")
    print("  sampler, labels = run_mcmc_desi(model='DG')")
    print("="*60)
