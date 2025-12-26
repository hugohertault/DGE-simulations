
! ============================================================
! DARK GEOMETRY PATCH FOR ECOSMOG/RAMSES
! ============================================================
! File: pm/poisson_commons.f90
! Add these routines

module dark_geometry
  use amr_parameters
  implicit none
  
  ! DG parameters (derived from first principles)
  real(dp), parameter :: dg_alpha_star = 0.075d0
  real(dp), parameter :: dg_beta = 0.666666666666667d0
  real(dp), parameter :: dg_k_j_eq = 0.05d0  ! h/Mpc
  real(dp), parameter :: dg_a_eq = 2.93255d-4
  
contains

  !------------------------------------------------------------
  ! Compute k_J(a) - Jeans scale
  !------------------------------------------------------------
  function dg_k_jeans(a) result(k_j)
    real(dp), intent(in) :: a
    real(dp) :: k_j
    real(dp) :: omega_m, omega_de, rho_ratio
    real(dp) :: factor_a, factor_eq
    
    omega_m = 0.3138d0
    omega_de = 0.6862d0
    
    rho_ratio = omega_m / omega_de * a**(-3)
    
    factor_a = rho_ratio**dg_beta - 1.0d0
    factor_eq = (omega_m/omega_de * dg_a_eq**(-3))**dg_beta - 1.0d0
    
    if (factor_a <= 0.0d0) then
      k_j = 0.001d0  ! DE regime
    else
      k_j = dg_k_j_eq * (a / dg_a_eq) * sqrt(factor_a / factor_eq)
    endif
    
    k_j = max(k_j, 0.001d0)
    k_j = min(k_j, 100.0d0)
    
  end function dg_k_jeans
  
  !------------------------------------------------------------
  ! Compute G_eff/G ratio
  !------------------------------------------------------------
  function dg_g_eff_ratio(k, a) result(g_ratio)
    real(dp), intent(in) :: k, a
    real(dp) :: g_ratio
    real(dp) :: k_j, k_ratio_sq
    
    k_j = dg_k_jeans(a)
    k_ratio_sq = (k / k_j)**2
    
    g_ratio = 1.0d0 + 2.0d0 * dg_alpha_star**2 / (1.0d0 + k_ratio_sq)
    
  end function dg_g_eff_ratio

end module dark_geometry

! ============================================================
! Modify Poisson solver (pm/poisson.f90)
! In the FFT-based solver, multiply by G_eff/G in k-space
! ============================================================

! After computing phi_k from delta_k:
!   phi_k = -4*pi*G * rho_bar * a^2 * delta_k / k^2
! 
! Apply DG modification:
!   phi_k = phi_k * dg_g_eff_ratio(k, aexp)

