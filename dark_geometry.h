/** @file dark_geometry.h 
 *  Header for Dark Geometry modified gravity implementation
 *  
 *  Dark Geometry (DG) unifies dark matter and dark energy as manifestations
 *  of a single scalar field - the Dark Boson - identified with the conformal
 *  mode of the spacetime metric.
 *  
 *  Key equation:
 *    G_eff(k,a)/G = 1 + 2*alpha_star(k)^2 / [1 + (k/k_J(a))^2]
 *  
 *  Parameters (derived from first principles):
 *    alpha_star_UV = 0.075 (Asymptotic Safety)
 *    beta = 2/3 (holographic area law)
 *    xi = 0.10 (non-minimal coupling for DG-E extension)
 *  
 *  Reference: Dark Geometry Mathematical Summary, December 2025
 */

#ifndef __DARK_GEOMETRY__
#define __DARK_GEOMETRY__

#include "common.h"

/* ========================================================================
 * FUNDAMENTAL DG PARAMETERS (derived, not free)
 * ======================================================================== */

/* UV coupling from Asymptotic Safety: alpha_star = g_star/(4 pi) sqrt(4/3) with g_star = 0.816 */
#define _DG_ALPHA_STAR_UV_ 0.075

/* IR coupling (after running) */
#define _DG_ALPHA_STAR_IR_ 0.083

/* Running coefficient */
#define _DG_BETA_ALPHA_ 0.01

/* Holographic exponent: beta = (d-1)/d = 2/3 for d=3 spatial dimensions */
#define _DG_BETA_ 0.666666666666667

/* Non-minimal coupling for DG-E: xi = beta/[4*(1+beta)] */
#define _DG_XI_ 0.10

/* Jeans scale at matter-radiation equality [h/Mpc] */
#define _DG_K_J_EQ_ 0.05

/* Scale factor at matter-radiation equality */
#define _DG_A_EQ_ 0.000293255

/* UV normalization scale [h/Mpc] */
#define _DG_K_UV_ 1.0

/* ========================================================================
 * DG PARAMETER STRUCTURE
 * ======================================================================== */

/**
 * Structure containing all Dark Geometry parameters
 */
struct dark_geometry_params {
  
  /* DG activation flags */
  short has_dg;              /* is Dark Geometry active? */
  short has_dg_e;            /* is DG-E extension (non-minimal coupling) active? */
  short has_dg_running;      /* is alpha_star running enabled? */
  
  /* Fundamental parameters (can be overridden for testing) */
  double alpha_star_UV;      /* UV coupling (default: 0.075) */
  double alpha_star_IR;      /* IR coupling (default: 0.083) */
  double beta_alpha;         /* running coefficient (default: 0.01) */
  double beta;               /* holographic exponent (default: 2/3) */
  double xi;                 /* non-minimal coupling (default: 0.10) */
  
  /* Scale parameters */
  double k_J_eq;             /* Jeans scale at equality [h/Mpc] */
  double a_eq;               /* scale factor at equality */
  double k_UV;               /* UV normalization scale [h/Mpc] */
  
  /* Cosmological parameters (copied from background) */
  double Omega_m_0;          /* matter density today */
  double Omega_DE_0;         /* dark energy density today */
  double h;                  /* reduced Hubble parameter */

};

/* ========================================================================
 * FUNCTION DECLARATIONS
 * ======================================================================== */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialize Dark Geometry parameters to default values
 */
int dark_geometry_init(struct dark_geometry_params * pdg);

/**
 * Free Dark Geometry structure
 */
int dark_geometry_free(struct dark_geometry_params * pdg);

/**
 * Compute alpha_star(k) with UV->IR running
 */
double dg_alpha_star(double k, struct dark_geometry_params * pdg);

/**
 * Compute Jeans scale k_J(a)
 */
double dg_k_Jeans(double a, struct dark_geometry_params * pdg);

/**
 * Compute effective mass squared m²_eff(a)
 */
double dg_mass_squared(double a, struct dark_geometry_params * pdg);

/**
 * Compute G_eff(k,a)/G - the key DG modification
 */
double dg_G_eff_ratio(double k, double a, struct dark_geometry_params * pdg);

/**
 * Compute power spectrum suppression factor S(k)
 */
double dg_suppression_factor(double k, double a, struct dark_geometry_params * pdg);

/**
 * Compute mu(k,a) for modified Poisson equation
 */
double dg_mu(double k, double a, struct dark_geometry_params * pdg);

/**
 * Compute eta(k,a) - gravitational slip
 */
double dg_eta(double k, double a, struct dark_geometry_params * pdg);

/**
 * Compute Sigma(k,a) - lensing modification
 */
double dg_Sigma(double k, double a, struct dark_geometry_params * pdg);

/**
 * Compute H(z) modification from DG-E extension
 */
double dg_H_modification(double z, struct dark_geometry_params * pdg);

/**
 * Check if DG modifications are active
 */
short dg_is_active(struct dark_geometry_params * pdg);

#ifdef __cplusplus
}
#endif

#endif /* __DARK_GEOMETRY__ */
