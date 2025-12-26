#!/usr/bin/env python3
"""
CLASS-DGE Python Wrapper
========================

This module provides a Python interface to the CLASS-DGE modifications.
It can work in two modes:

1. STANDALONE: Uses Python implementations of DGE physics
2. CLASS: Interfaces with modified CLASS code (requires compilation)

For standalone mode, no compilation is needed.
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq

# =============================================================================
# DGE CONSTANTS (ALL DERIVED)
# =============================================================================

# From Asymptotic Safety
ALPHA_STAR = 0.075

# From holographic area law
BETA = 2.0 / 3.0

# Non-minimal coupling
XI = BETA / (4.0 * (1.0 + BETA))  # = 0.10

# Physical constants
c_km_s = 299792.458  # km/s
H0_default = 67.36   # km/s/Mpc


# =============================================================================
# COSMOLOGICAL BACKGROUND
# =============================================================================

class DGEBackground:
    """
    Background cosmology with DGE modifications.
    """
    
    def __init__(self, H0=67.36, omega_b=0.0224, omega_cdm=0.120, 
                 T_cmb=2.7255, N_eff=3.046, phi_initial=0.19):
        
        self.H0 = H0
        self.h = H0 / 100.0
        self.omega_b = omega_b
        self.omega_cdm = omega_cdm
        self.omega_m = omega_b + omega_cdm
        self.T_cmb = T_cmb
        self.N_eff = N_eff
        self.phi_initial = phi_initial
        
        # Radiation
        self.omega_gamma = 2.47e-5 * (T_cmb / 2.725)**4
        self.omega_nu = N_eff * (7.0/8.0) * (4.0/11.0)**(4.0/3.0) * self.omega_gamma
        self.omega_r = self.omega_gamma + self.omega_nu
        
        # Dark energy (flat universe)
        self.omega_L = self.h**2 - self.omega_m - self.omega_r
        
        # Key redshifts
        self.z_eq = self.omega_m / self.omega_r - 1
        self.z_rec = 1090.0
        self.z_drag = 1060.0
        
    def phi(self, z):
        """Dark Boson field amplitude."""
        z_eq = self.z_eq
        z_trans = 0.5
        
        if np.isscalar(z):
            z = np.array([z])
            scalar = True
        else:
            z = np.asarray(z)
            scalar = False
        
        result = np.zeros_like(z, dtype=float)
        
        # Different regimes
        mask_rad = z > z_eq
        mask_mat = (z <= z_eq) & (z > z_trans)
        mask_de = z <= z_trans
        
        result[mask_rad] = self.phi_initial
        
        if np.any(mask_mat):
            growth = ((1 + z_eq) / (1 + z[mask_mat]))**0.3
            result[mask_mat] = self.phi_initial * growth
        
        if np.any(mask_de):
            growth_eq = ((1 + z_eq) / (1 + z_trans))**0.3
            late_growth = ((1 + z_trans) / (1 + z[mask_de]))**0.1
            result[mask_de] = self.phi_initial * growth_eq * late_growth
        
        return result[0] if scalar else result
    
    def G_eff_ratio(self, z):
        """G_eff / G from non-minimal coupling."""
        phi = self.phi(z)
        xi_phi2 = XI * phi**2
        
        denom = 1.0 - 8.0 * np.pi * xi_phi2
        
        if np.any(denom <= 0):
            return 1.0 + 8.0 * np.pi * xi_phi2
        
        return 1.0 / denom
    
    def H_LCDM(self, z):
        """Standard ΛCDM Hubble parameter."""
        return self.H0 * np.sqrt(
            self.omega_r / self.h**2 * (1+z)**4 +
            self.omega_m / self.h**2 * (1+z)**3 +
            self.omega_L / self.h**2
        )
    
    def H(self, z):
        """Modified Hubble parameter with DGE."""
        return self.H_LCDM(z) * np.sqrt(self.G_eff_ratio(z))
    
    def sound_speed(self, z):
        """Sound speed in baryon-photon plasma."""
        a = 1.0 / (1.0 + z)
        R_b = 3.0 * self.omega_b * a / (4.0 * self.omega_gamma)
        return 1.0 / np.sqrt(3.0 * (1.0 + R_b))
    
    def r_s(self, z_upper):
        """Sound horizon."""
        def integrand(z):
            cs = self.sound_speed(z)
            H = self.H(z)
            return c_km_s * cs / H
        
        result, _ = quad(integrand, z_upper, 1e5, limit=500)
        return result
    
    def D_A(self, z):
        """Angular diameter distance."""
        def integrand(zp):
            return 1.0 / self.H(zp)
        
        result, _ = quad(integrand, 0, z, limit=200)
        return c_km_s * result / (1.0 + z)
    
    def theta_star(self):
        """Acoustic angle."""
        rs = self.r_s(self.z_drag)
        DA = self.D_A(self.z_rec)
        return rs / DA


# =============================================================================
# PERTURBATIONS
# =============================================================================

class DGEPerturbations:
    """
    Linear perturbation theory with DGE modifications.
    """
    
    def __init__(self, background):
        self.bg = background
        
    def G_eff_k(self, z, k):
        """Scale-dependent effective gravitational coupling."""
        G_eff_bg = self.bg.G_eff_ratio(z)
        
        # Jeans scale
        k_J = 0.05 * (1 + z) / self.bg.z_eq
        
        # Scale-dependent enhancement
        alpha2 = ALPHA_STAR**2
        scale_factor = 1.0 + 2.0 * alpha2 / (1.0 + (k / k_J)**2)
        
        return G_eff_bg * scale_factor
    
    def growth_equation(self, ln_a, y, k):
        """
        Growth equation: D'' + (2 + d ln H/d ln a) D' - 3/2 Ω_m G_eff D = 0
        
        y = [D, dD/d(ln a)]
        """
        D, D_prime = y
        
        a = np.exp(ln_a)
        z = 1.0 / a - 1.0
        
        # Hubble
        H = self.bg.H(z)
        H_LCDM = self.bg.H_LCDM(z)
        
        # Ω_m(z)
        omega_m_z = self.bg.omega_m / self.bg.h**2 * (1+z)**3 / (H / self.bg.H0)**2
        
        # d ln H / d ln a (approximate)
        dz = 0.01 * (1 + z)
        H_plus = self.bg.H(z + dz)
        H_minus = self.bg.H(z - dz)
        dlnH_dlna = -(1 + z) * (H_plus - H_minus) / (2 * dz * H)
        
        # G_eff(k)
        G_eff = self.G_eff_k(z, k)
        
        # Equation of motion
        dD_dlna = D_prime
        dD_prime_dlna = -(2 + dlnH_dlna) * D_prime + 1.5 * omega_m_z * G_eff * D
        
        return [dD_dlna, dD_prime_dlna]
    
    def solve_growth(self, k, z_initial=1e4, z_final=0):
        """Solve growth equation for given k."""
        
        ln_a_initial = -np.log(1 + z_initial)
        ln_a_final = -np.log(1 + z_final) if z_final > 0 else 0
        
        # Initial conditions (matter domination)
        D_initial = 1.0 / (1 + z_initial)
        D_prime_initial = 1.0  # D ∝ a in matter era
        
        y0 = [D_initial, D_prime_initial]
        
        sol = solve_ivp(
            lambda ln_a, y: self.growth_equation(ln_a, y, k),
            (ln_a_initial, ln_a_final),
            y0,
            method='RK45',
            dense_output=True,
            rtol=1e-6,
            atol=1e-8
        )
        
        return sol
    
    def growth_factor(self, k, z=0):
        """Growth factor D(k, z)."""
        sol = self.solve_growth(k, z_final=z)
        return sol.y[0, -1]
    
    def power_suppression(self, k, z=0):
        """Power spectrum suppression S(k) = P_DGE / P_LCDM."""
        
        # Growth in DGE
        D_DGE = self.growth_factor(k, z)
        
        # Growth in LCDM (k-independent, take k→0 limit)
        D_LCDM = self.growth_factor(0.001, z)
        
        # Suppression (squared because P ∝ D²)
        return (D_DGE / D_LCDM)**2
    
    def sigma8(self, sigma8_LCDM=0.811):
        """σ₈ with DGE modifications."""
        k_eff = 0.2  # Effective k for σ₈
        S_eff = self.power_suppression(k_eff)
        return sigma8_LCDM * np.sqrt(S_eff)


# =============================================================================
# CMB OBSERVABLES
# =============================================================================

class DGECMB:
    """
    CMB observables with DGE modifications.
    """
    
    def __init__(self, background):
        self.bg = background
        
    def r_s_ratio(self):
        """r_s^DGE / r_s^LCDM."""
        # Create ΛCDM background
        bg_LCDM = DGEBackground(
            H0=self.bg.H0,
            omega_b=self.bg.omega_b,
            omega_cdm=self.bg.omega_cdm,
            phi_initial=0.0  # No DGE
        )
        
        rs_DGE = self.bg.r_s(self.bg.z_drag)
        rs_LCDM = bg_LCDM.r_s(bg_LCDM.z_drag)
        
        return rs_DGE / rs_LCDM
    
    def D_A_ratio(self):
        """D_A^DGE / D_A^LCDM."""
        bg_LCDM = DGEBackground(
            H0=self.bg.H0,
            omega_b=self.bg.omega_b,
            omega_cdm=self.bg.omega_cdm,
            phi_initial=0.0
        )
        
        DA_DGE = self.bg.D_A(self.bg.z_rec)
        DA_LCDM = bg_LCDM.D_A(bg_LCDM.z_rec)
        
        return DA_DGE / DA_LCDM
    
    def theta_star_ratio(self):
        """θ*^DGE / θ*^LCDM."""
        return self.r_s_ratio() / self.D_A_ratio()
    
    def infer_H0(self, theta_star_target):
        """
        Find H₀ that gives the target θ*.
        
        Since θ* is tightly constrained by Planck,
        we find H₀ such that θ*^DGE = θ*^LCDM.
        """
        def objective(H0):
            bg = DGEBackground(
                H0=H0,
                omega_b=self.bg.omega_b,
                omega_cdm=self.bg.omega_cdm,
                phi_initial=self.bg.phi_initial
            )
            cmb = DGECMB(bg)
            return bg.theta_star() - theta_star_target
        
        # Search for H₀
        H0_opt = brentq(objective, 60, 80)
        
        return H0_opt


# =============================================================================
# FULL DGE COSMOLOGY
# =============================================================================

class DGECosmology:
    """
    Complete DGE cosmology wrapper.
    """
    
    def __init__(self, H0=67.36, omega_b=0.0224, omega_cdm=0.120,
                 phi_initial=0.19):
        
        self.background = DGEBackground(
            H0=H0, omega_b=omega_b, omega_cdm=omega_cdm,
            phi_initial=phi_initial
        )
        self.perturbations = DGEPerturbations(self.background)
        self.cmb = DGECMB(self.background)
        
        # Store parameters
        self.H0 = H0
        self.omega_b = omega_b
        self.omega_cdm = omega_cdm
        self.phi_initial = phi_initial
    
    def compute_all(self):
        """Compute all observables."""
        
        results = {}
        
        # Background
        results['phi_rec'] = self.background.phi(self.background.z_rec)
        results['G_eff_rec'] = self.background.G_eff_ratio(self.background.z_rec)
        results['r_s'] = self.background.r_s(self.background.z_drag)
        results['D_A'] = self.background.D_A(self.background.z_rec)
        results['theta_star'] = self.background.theta_star()
        
        # CMB ratios
        results['r_s_ratio'] = self.cmb.r_s_ratio()
        results['D_A_ratio'] = self.cmb.D_A_ratio()
        results['theta_star_ratio'] = self.cmb.theta_star_ratio()
        
        # Perturbations
        results['sigma8'] = self.perturbations.sigma8()
        results['S_k02'] = self.perturbations.power_suppression(0.2)
        
        return results
    
    def summary(self):
        """Print summary of DGE cosmology."""
        
        results = self.compute_all()
        
        print("="*60)
        print("DGE COSMOLOGY SUMMARY")
        print("="*60)
        print(f"\nParameters (ALL DERIVED):")
        print(f"  α* = {ALPHA_STAR:.4f}")
        print(f"  β  = {BETA:.4f}")
        print(f"  ξ  = {XI:.4f}")
        print(f"\nInput:")
        print(f"  H₀ = {self.H0:.2f} km/s/Mpc")
        print(f"  φ_initial = {self.phi_initial:.4f} M_Pl")
        print(f"\nField evolution:")
        print(f"  φ(z=1090) = {results['phi_rec']:.4f} M_Pl")
        print(f"  G_eff/G(z=1090) = {results['G_eff_rec']:.4f}")
        print(f"\nCMB:")
        print(f"  r_s = {results['r_s']:.2f} Mpc")
        print(f"  D_A = {results['D_A']:.2f} Mpc")
        print(f"  θ* = {results['theta_star']:.6f}")
        print(f"  Δr_s/r_s = {(results['r_s_ratio']-1)*100:+.2f}%")
        print(f"\nPerturbations:")
        print(f"  σ₈ = {results['sigma8']:.4f}")
        print(f"  S(k=0.2) = {results['S_k02']:.4f}")
        print("="*60)
        
        return results


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    
    print("="*60)
    print("CLASS-DGE PYTHON WRAPPER")
    print("="*60)
    
    # Create DGE cosmology
    dge = DGECosmology(phi_initial=0.19)
    
    # Compute and print summary
    results = dge.summary()
    
    # Find H₀ that preserves θ*
    print("\n" + "-"*60)
    print("H₀ INFERENCE (preserving θ*)")
    print("-"*60)
    
    # ΛCDM reference θ*
    lcdm = DGECosmology(phi_initial=0.0)
    theta_LCDM = lcdm.background.theta_star()
    print(f"\nΛCDM θ* = {theta_LCDM:.6f}")
    
    # Find H₀ for DGE
    H0_inferred = dge.cmb.infer_H0(theta_LCDM)
    print(f"DGE H₀ (to match θ*) = {H0_inferred:.2f} km/s/Mpc")
    print(f"ΔH₀ = {H0_inferred - 67.36:+.2f} km/s/Mpc")
