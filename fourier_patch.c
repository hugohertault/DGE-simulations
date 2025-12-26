/**
 * Dark Geometry Patch for CLASS fourier.c
 * ========================================
 * 
 * This file contains the modifications to apply to CLASS's fourier.c
 * to implement the Dark Geometry power spectrum suppression.
 * 
 * Location: Around line 1380 in fourier.c, after P(k) is computed
 * 
 * INSTRUCTIONS:
 * 1. Open class_public/source/fourier.c
 * 2. Find the section where pk_nl is computed (search for "nonlinear_method")
 * 3. Add the DG suppression AFTER the standard P(k) calculation
 * 4. Recompile CLASS with: make clean && make
 */

/* =========================================================
 * ADD THIS CODE BLOCK after pk_nl is computed
 * ========================================================= */

/* Dark Geometry: Apply S(k) suppression to power spectrum */
#ifdef _DARK_GEOMETRY_

#include "dark_geometry.h"

/* Apply DG suppression after computing P(k) */
{
    double k_J, S_k, x_sq;
    double a = 1.0 / (1.0 + z);  /* scale factor */
    
    /* Compute Jeans scale at this epoch */
    k_J = dg_k_jeans(a, pba->Omega0_cdm + pba->Omega0_b, 
                     1.0 - pba->Omega0_cdm - pba->Omega0_b);
    
    /* Compute suppression factor */
    x_sq = (k / k_J) * (k / k_J);
    S_k = (1.0 + DG_S_MAX * x_sq) / (1.0 + x_sq);
    
    /* Apply suppression to both linear and non-linear P(k) */
    pk_l *= S_k;
    pk_nl *= S_k;
    
    /* Debug output (comment out for production) */
    /* printf("DG: k=%.4f, k_J=%.4f, S=%.4f\n", k, k_J, S_k); */
}

#endif /* _DARK_GEOMETRY_ */


/* =========================================================
 * ALTERNATIVE: Modify the perturbations via G_eff
 * This approach modifies the growth rate directly
 * ========================================================= */

/* In perturbations.c, modify the Poisson equation:
 * 
 * Original:
 *   k^2 phi = -4*pi*G * rho * delta
 * 
 * With DG:
 *   k^2 phi = -4*pi*G_eff(k,a) * rho * delta
 * 
 * where G_eff/G = 1 + 2*alpha_star^2 / (1 + (k/k_J)^2)
 */

/* =========================================================
 * COMPILATION FLAGS
 * ========================================================= */

/* Add to Makefile:
 * 
 * CFLAGS += -D_DARK_GEOMETRY_
 * 
 * Or compile with:
 * make CFLAGS="-D_DARK_GEOMETRY_"
 */


/* =========================================================
 * FULL FUNCTION: dg_apply_suppression
 * Call this after primordial_spectrum_at_k() returns
 * ========================================================= */

/**
 * Apply Dark Geometry suppression to power spectrum
 * 
 * @param k          Wavenumber in h/Mpc
 * @param a          Scale factor
 * @param pk_in      Input power spectrum
 * @param Omega_m    Matter density parameter
 * @return           Suppressed power spectrum
 */
double dg_apply_suppression(double k, double a, double pk_in, double Omega_m) {
    
    /* DG parameters (derived from first principles) */
    const double alpha_star = 0.075;
    const double S_max = 0.882;
    const double k_J_eq = 0.05;  /* h/Mpc at equality */
    const double a_eq = 2.93e-4;
    const double beta = 2.0/3.0;
    
    double Omega_de = 1.0 - Omega_m;
    double rho_ratio, k_J, x_sq, S_k;
    
    /* Compute k_J(a) */
    rho_ratio = Omega_m / Omega_de * pow(a, -3);
    
    if (rho_ratio <= 1.0) {
        /* DE dominated: minimal suppression */
        k_J = 0.001;
    } else {
        double factor_a = pow(rho_ratio, beta) - 1.0;
        double factor_eq = pow(Omega_m/Omega_de * pow(a_eq, -3), beta) - 1.0;
        k_J = k_J_eq * (a / a_eq) * sqrt(factor_a / factor_eq);
    }
    
    /* Ensure k_J is bounded */
    if (k_J < 0.001) k_J = 0.001;
    if (k_J > 100.0) k_J = 100.0;
    
    /* Compute S(k) */
    x_sq = (k / k_J) * (k / k_J);
    S_k = (1.0 + S_max * x_sq) / (1.0 + x_sq);
    
    return pk_in * S_k;
}
