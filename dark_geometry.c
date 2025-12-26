/** @file dark_geometry.c 
 *  Implementation of Dark Geometry modified gravity functions
 *  
 *  Key insight from the theory:
 *  - G_eff/G = 1 + 2*alpha^2/(1+(k/k_J)^2) is ALWAYS >= 1
 *  - The SUPPRESSION comes from differential growth over N~8 e-folds
 *  - At k << k_J: faster growth (more boost)
 *  - At k >> k_J: slower growth (less boost, approaches GR)
 *  - Net effect: P(k>>k_J)/P(k<<k_J) ~ 0.88 (12% suppression)
 *  
 *  Reference: Dark Geometry Mathematical Summary, December 2025
 */

#include "dark_geometry.h"
#include <math.h>

/* ============================================================================
 * CONSTANTS derived from theory
 * ============================================================================ */

/* Number of e-folds from equality to today: N = ln(a_0/a_eq) = ln(3409) */
static const double N_EFOLDS = 8.13;

/* Maximum suppression factor: S_max = exp(-2*N*Delta_n) = 0.882 */
static const double S_MAX = 0.882;

/* ============================================================================
 * INITIALIZATION
 * ============================================================================ */

int dark_geometry_init(struct dark_geometry_params * pdg) {
  
  /* Default: DG active */
  pdg->has_dg = _TRUE_;
  pdg->has_dg_e = _TRUE_;
  pdg->has_dg_running = _TRUE_;
  
  /* Fundamental parameters (derived from first principles) */
  pdg->alpha_star_UV = _DG_ALPHA_STAR_UV_;
  pdg->alpha_star_IR = _DG_ALPHA_STAR_IR_;
  pdg->beta_alpha = _DG_BETA_ALPHA_;
  pdg->beta = _DG_BETA_;
  pdg->xi = _DG_XI_;
  
  /* Scale parameters */
  pdg->k_J_eq = _DG_K_J_EQ_;
  pdg->a_eq = _DG_A_EQ_;
  pdg->k_UV = _DG_K_UV_;
  
  /* Cosmological parameters (will be set from background) */
  pdg->Omega_m_0 = 0.315;
  pdg->Omega_DE_0 = 0.685;
  pdg->h = 0.6736;
  
  return _SUCCESS_;
}

int dark_geometry_free(struct dark_geometry_params * pdg) {
  return _SUCCESS_;
}

/* ============================================================================
 * CORE DG FUNCTIONS
 * ============================================================================ */

/**
 * Compute alpha_star(k) with UV->IR running
 */
double dg_alpha_star(double k, struct dark_geometry_params * pdg) {
  
  double alpha;
  
  if (pdg->has_dg == _FALSE_) {
    return 0.0;
  }
  
  if (pdg->has_dg_running == _FALSE_) {
    return pdg->alpha_star_UV;
  }
  
  if (k <= 1.e-10) {
    k = 1.e-10;
  }
  
  alpha = pdg->alpha_star_UV * (1.0 + pdg->beta_alpha * log(pdg->k_UV / k));
  
  if (alpha < 0.0) alpha = 0.0;
  if (alpha > 0.15) alpha = 0.15;
  
  return alpha;
}

/**
 * Compute Jeans scale k_J(a)
 */
double dg_k_Jeans(double a, struct dark_geometry_params * pdg) {
  
  double k_J;
  double rho_ratio_0, rho_ratio_a, rho_ratio_eq;
  double factor_a, factor_eq;
  
  if (pdg->has_dg == _FALSE_) {
    return 1.e10;
  }
  
  rho_ratio_0 = pdg->Omega_m_0 / pdg->Omega_DE_0;
  rho_ratio_a = rho_ratio_0 * pow(a, -3.0);
  rho_ratio_eq = rho_ratio_0 * pow(pdg->a_eq, -3.0);
  
  factor_a = pow(rho_ratio_a, pdg->beta) - 1.0;
  factor_eq = pow(rho_ratio_eq, pdg->beta) - 1.0;
  
  /* Handle DE regime */
  if (factor_a <= 0.0) {
    double a_trans = pow(rho_ratio_0, 1.0/3.0);
    double rho_ratio_trans = 1.0;
    double factor_trans = pow(rho_ratio_trans, pdg->beta) - 1.0;
    
    if (factor_trans <= 0.01) factor_trans = 0.01;
    if (factor_eq <= 0.01) factor_eq = 0.01;
    
    k_J = pdg->k_J_eq * (a_trans / pdg->a_eq) * sqrt(factor_trans / factor_eq);
    return (k_J > 0.001) ? k_J : 0.001;
  }
  
  if (factor_eq <= 0.0) {
    return pdg->k_J_eq;
  }
  
  k_J = pdg->k_J_eq * (a / pdg->a_eq) * sqrt(factor_a / factor_eq);
  
  if (k_J < 0.001) k_J = 0.001;
  if (k_J > 100.0) k_J = 100.0;
  
  return k_J;
}

/**
 * Compute G_eff(k,a)/G - ALWAYS >= 1
 */
double dg_G_eff_ratio(double k, double a, struct dark_geometry_params * pdg) {
  
  double alpha, k_J;
  double k_ratio_sq;
  double G_ratio;
  
  if (pdg->has_dg == _FALSE_) {
    return 1.0;
  }
  
  alpha = dg_alpha_star(k, pdg);
  k_J = dg_k_Jeans(a, pdg);
  
  k_ratio_sq = (k / k_J) * (k / k_J);
  
  G_ratio = 1.0 + 2.0 * alpha * alpha / (1.0 + k_ratio_sq);
  
  return G_ratio;
}

/**
 * Compute power spectrum suppression S(k)
 * 
 * S(k) = 1 - (1 - S_max) * [1 - 1/(1 + (k/k_J)^2)]
 *      = (1 + S_max * x^2) / (1 + x^2)  where x = k/k_J
 */
double dg_suppression_factor(double k, double a, struct dark_geometry_params * pdg) {
  
  double k_J, x2;
  double S;
  
  if (pdg->has_dg == _FALSE_) {
    return 1.0;
  }
  
  k_J = dg_k_Jeans(a, pdg);
  x2 = (k / k_J) * (k / k_J);
  
  /* S = 1 - (1-S_MAX) * x2/(1+x2) = (1 + S_MAX*x2) / (1+x2) */
  S = (1.0 + S_MAX * x2) / (1.0 + x2);
  
  return S;
}

/**
 * Compute mu(k,a) - effective modification for Poisson equation
 * 
 * NOTE: For the robust implementation, suppression is applied directly
 * to P(k) in fourier.c, so mu = 1 here to avoid double counting.
 */
double dg_mu(double k, double a, struct dark_geometry_params * pdg) {
  
  /* Suppression applied in fourier.c, so mu = 1 here */
  return 1.0;
}

/**
 * eta = Psi/Phi = 1 (no slip in DG)
 */
double dg_eta(double k, double a, struct dark_geometry_params * pdg) {
  return 1.0;
}

/**
 * Sigma for lensing
 */
double dg_Sigma(double k, double a, struct dark_geometry_params * pdg) {
  double mu = dg_mu(k, a, pdg);
  double eta = dg_eta(k, a, pdg);
  
  return mu * (1.0 + eta) / 2.0;
}

/**
 * H(z) modification from DG-E
 */
double dg_H_modification(double z, struct dark_geometry_params * pdg) {
  
  double H_ratio;
  double eta_H = 0.8;
  double z_rec = 1090.0;
  double sigma_z = 500.0;
  double weight;
  
  if (pdg->has_dg_e == _FALSE_) {
    return 1.0;
  }
  
  weight = exp(-(z - z_rec) * (z - z_rec) / (2.0 * sigma_z * sigma_z));
  H_ratio = 1.0 + eta_H * pdg->xi * weight;
  
  return H_ratio;
}

/**
 * Check if DG is active
 */
short dg_is_active(struct dark_geometry_params * pdg) {
  if (pdg->has_dg == _TRUE_ && pdg->alpha_star_UV > 0.0) {
    return _TRUE_;
  }
  return _FALSE_;
}
