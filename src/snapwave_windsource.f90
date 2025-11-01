module snapwave_windsource
!
implicit none
!
contains
!
subroutine numerical_limiter(ntheta, ee, aa, waveps, hh, dtheta, rho, g, gamma, sigmin, sigmax, H, E, A, sig)
   implicit none
   !
   integer, intent(in) :: ntheta
   real*4, dimension(ntheta), intent(inout) :: ee
   real*4, dimension(ntheta), intent(inout) :: aa
   real*4, intent(in) :: waveps
   real*4, intent(in) :: hh
   real*4, intent(in) :: dtheta
   real*4, intent(in) :: rho
   real*4, intent(in) :: g
   real*4, intent(in) :: gamma
   real*4, intent(in) :: sigmin
   real*4, intent(in) :: sigmax
   !
   real*4, intent(out) :: H
   real*4, intent(out) :: E
   real*4, intent(out) :: A
   real*4, intent(out) :: sig

   integer itheta
   real*4 :: depthlimfac
   real*4, dimension(:), allocatable :: sigt
   !
   allocate (sigt(ntheta))
   !
   ee(:) = max(ee(:), waveps)
   aa(:) = max(ee(:), waveps) / sigmax !sigmax
   !
   sigt = ee / aa
   E = sum(ee(:)) * dtheta
   A = sum(aa(:)) * dtheta
   !
   ! depth limitation of energy and corresponding wave period, from directionally integrated params
   !
   H = sqrt(8.0 * E / rho / g)
   depthlimfac = max(1.0, (H / (gamma * hh))**2)
   H = min(H, gamma * hh)
   !
   E = E / depthlimfac
   A = A / depthlimfac
   !
   do itheta = 1, ntheta
      ee(itheta) = ee(itheta) / depthlimfac
      aa(itheta) = aa(itheta) / depthlimfac
   end do
   !
   ! limit frequency to range
   !
   sig = E / A
   sig = max(sig, sigmin)
   sig = min(sig, sigmax)
   do itheta = 1, ntheta
      aa(itheta) = ee(itheta) / sig
   end do
   A = sum(aa(:)) * dtheta
   sig = E / A
   !
end subroutine numerical_limiter
!
subroutine compute_celerities(hh, sig, sinth, costh, ntheta, gamma, dhdx, dhdy, sinhkh, Hmx, kwav, cg, ctheta)
   use snapwave_data, only: pi
   implicit none
   integer, intent(in) :: ntheta
   real*4, intent(in) :: hh, sig, dhdx, dhdy
   real*4, dimension(ntheta), intent(in) :: sinth, costh
   real*4, intent(in) :: gamma
   !
   real*4, intent(out) :: sinhkh, Hmx, cg, kwav
   real*4, dimension(ntheta), intent(out) :: ctheta
   !
   integer :: itheta
   !
   ! Compute celerities and refraction speed
   ! in a seperate routine because will be called more often in a model with varying wave period
   call disper_nr(hh, 2.0 * pi / sig, kwav, cg) ! 15% faster
   !
   sinhkh = sinh(min(kwav * hh, 50.0))
   !Hmx =0.88/kwav*tanh(gamma*kwav*hh/0.88)
   Hmx = gamma * hh
   !
   do itheta = 1, ntheta
      ctheta(itheta) = sig / sinh(min(2.0 * kwav * hh, 50.0)) * (dhdx * sinth(itheta) - dhdy * costh(itheta))
      ! Limit unrealistic refraction speed to 1/2 pi per wave period
      ctheta(itheta) = sign(1.0, ctheta(itheta)) * min(abs(ctheta(itheta)), sig / 4.0)
   end do
   !
end subroutine compute_celerities
!
subroutine windinput(u10, rho, g, hh, ntheta, windspreadfac, E, A, cg, eeprev, aaprev, ds, WsorE, WsorA, jadcgdx)
   implicit none
   real*4, intent(in) :: u10 !< [m/s] wind speed
   real*4, intent(in) :: rho !< [kg/m3] density of water
   real*4, intent(in) :: g !< [m/s2] gravitational acceleration
   real*4, intent(in) :: hh !< [m] water depth
   integer, intent(in) :: ntheta !< [-] no. directional bins
   real*4, dimension(ntheta), intent(in) :: windspreadfac !< [-] distribution array for wind input
   !
   real*4, intent(in) :: E !< [J/m2] nodal wave energy
   real*4, intent(in) :: A !< [s] nodal representative wave period
   real*4, intent(in) :: cg !< [s] nodal group velocity
   real*4, dimension(ntheta), intent(in) :: eeprev !< [J/s/m2] nodal wave action at upwind point for wind input
   real*4, dimension(ntheta), intent(in) :: aaprev !< [J/m2] nodal waveenergy at upwind point for wind input
   real*4, dimension(ntheta), intent(in) :: ds !< [m] distance to upwind point for wind input
   !
   real*4, dimension(ntheta), intent(out) :: wsorE !< [J/rad/s]  wind input energy per second
   real*4, dimension(ntheta), intent(out) :: wsorA !< [J/rad/s/s]  wind input action per second
   !
   integer, intent(in) :: jadcgdx
   !
   ! Local variables and arrays
   !
   real*4 :: dTdx !  [s/m] gradient of wave period
   real*4 :: T
   real*4 :: Tprev !  [s/m] auxiliary for gradient of T
   real*4 :: deltaT !  [s/m] auxiliary for gradient of T
   real*4 :: dcgdT !  [m/s/s] rate of change of group velocity with wave period

   real*4 :: Eful = 0.0036 !  [-] fully developed dimensionless wave energy (Pierson Moskowitz 1964)
   real*4 :: Tful = 7.69 !  [-] fully developed dimensionless peak period (Pierson Moskowitz 1964)
   real*4 :: aa1 = 0.00288 !  [-] shape parameter wave growth curves (Kahma Calkoen (1992))
   real*4 :: bb1 = 0.45 !  [-] shape parameter wave growth curves (Kahma Calkoen (1992))
   real*4 :: aa2 = 0.459 !  [-] shape parameter wave growth curves (Kahma Calkoen (1992))
   real*4 :: bb2 = 0.27 !  [-] shape parameter wave growth curves (Kahma Calkoen (1992))
   real*4 :: aa3 = 0.13 !  [-] shape parameter wave growth curves (Breugem Holthuijzen (2007))
   real*4 :: bb3 = 0.65 !  [-] shape parameter wave growth curves (Breugem Holthuijzen (2007))
   real*4 :: aa4 = 5.0 !  [-] shape parameter wave growth curves (Breugem Holthuijzen (2007))
   real*4 :: bb4 = 0.375 !  [-] shape parameter wave growth curves (Breugem Holthuijzen (2007))

   real*4 :: Edmlss, dEtmp
   real*4 :: Tdmlss
   real*4 :: cgdmlss
   real*4 :: ddmlss
   real*4 :: Emaxddmlss
   real*4 :: Tmaxddmlss
   real*4 :: wsorEdlss
   real*4 :: wsorAdlss

   real*4 :: pi, fE, fT, dE, dT, dcgdx
   real*4 :: cg1, cg2, kdum ! auxiliary variables to determine dcg/dT
   integer :: itheta
   !
   pi = 4.d0 * atan(1.d0)
   !
   !compute dimensionless wave state, maximized by Breugem and Holthuizen fully developed sea states
   !

   ! as a limiter on the wave input, we maximize the energy and wave period with which the source term is computed.
   ! this means there is still wind input, but its amount is bounded by the water depth
   ddmlss = g * hh / u10**2
   Emaxddmlss = min(aa3**2 / 16.0 * ddmlss**(2 * bb3), Eful)
   Tmaxddmlss = min(aa4 * ddmlss**bb4, Tful) !2.0d0*aa2*(16.0d0*Emaxddmlss/aa1**2)**(bb2/2.0d0/bb1) ! this is twice the wave period that would correspond to the maximum depth limited energy
   !
   cgdmlss = abs(cg) / u10 !0.8*abs(cg) / u10
   Edmlss = min(E * g / rho / u10**4, Emaxddmlss)
   !T = min(2.0d0*pi*A / E, 2*aa2*(16*E/aa1**2)**(bb2/2/bb1) )
   T = 2.0 * pi * A / E
   Tdmlss = min(g * T / u10, Tmaxddmlss)
   !
   ! dimensionless magnitude of source terms, based on Kahma and Calkoen
   !
   fE = 16.d0 / (2.0 * bb1 * aa1 * aa1) * (16.0 * Edmlss / aa1 / aa1)**(0.5 / bb1 - 1.0)
   dE = cgdmlss / fE
   !
   fT = 1.d0 / aa2 / bb2 * (Tdmlss / aa2)**(1.0 / bb2 - 1)
   dT = cgdmlss / fT
   !
   do itheta = 1, ntheta
      !
      !gradT component computed from dimensional parameters, made dimensionless in last step
      !
      Tprev = 2.d0 * pi * aaprev(itheta) / max(eeprev(itheta),0.001)
      Tprev=max(Tprev,1.)
      deltaT = T - Tprev
      call disper_nr(hh, T, kdum, cg1)
      call disper_nr(hh, Tprev, kdum, cg2)
      if (abs(deltaT) > 1e-6) then
         dcgdT = (cg1 - cg2) / deltaT
      else
         dcgdT = 0.
      end if
      dTdx = jadcgdx * max((T - Tprev) / ds(itheta), 0.0)
      dcgdx = E * dcgdT * dTdx / u10**3 / rho !  dimensionless
      dcgdx = min(dE, abs(dcgdx))
      !
      !dimensionalize and distribute over directional bins
      !
      !dT = dT * max(0.0d0,tanh(2*(Tful-Tdmlss)/Tful)) !used to be 8.5 for 1D model. Here we make it 1000 to make sure it works with 36 dirs
      dT = dT * max(0.0, tanh(1000.0 * (Tful - Tdmlss) / Tful)) !used to be 8.5 for 1D model. Here we make it 1000 to make sure it works with 36 dirs
      !dT = min(dT, oneoverdt * Tmaxddmlss/2) !numerically limit the amount of input per timestep to maximally half the fully developed conditions
      !dT = min(dT, oneoverdt * Tdmlss/10) !numerically limit the amount of input per timestep to maximally half the fully developed conditions
      !dE = min(dE, oneoverdt * Emaxddmlss/2) !numerically limit the amount of input per timestep to maximally half the fully developed conditions
      !dE = min(dE, oneoverdt * Edmlss/10) !numerically limit the amount of input per timestep to maximally half the fully developed conditions

      dEtmp = (dE + dcgdx) * max(0.0, tanh(2.0 * (Eful - Edmlss) / Eful))
      wsorAdlss = windspreadfac(itheta) * 0.5 / pi * (Tdmlss * dEtmp + Edmlss * dT) !
      wsorEdlss = windspreadfac(itheta) * dEtmp
      !
      ! make dimensional growth rates
      !
      wsorE(itheta) = max(u10**3 * rho * wsorEdlss, 0.0) !
      wsorA(itheta) = max(u10**4 * rho / g * wsorAdlss, 0.0)
      !SwT = dT
   end do
   !
end subroutine windinput
!
subroutine disper_nr(h, T, k, cg)
   use snapwave_data, only: pi, g
   implicit none
   real*4, intent(in) :: h ! water depth (m)
   real*4, intent(in) :: T ! wave period (s)
   real*4, intent(out) :: k ! wavenumber (rad/m)
   real*4, intent(out) :: cg ! wavenumber (rad/m)

   real*4, parameter :: tol = 1.0d-5 ! convergence tolerance
   real*4 :: omega, kh, f, df, k_old, c, arg, n
   integer :: i, max_iter = 20

   ! Calculate angular frequency from period
   omega = 2.0d0 * pi / T

   ! Initial guess using deep water approximation: k = ω²/g
   k = omega * omega / g

   ! Newton-Raphson iteration
   do i = 1, max_iter
      kh = k * h
      f = g * k * tanh(kh) - omega * omega ! residual
      df = g * (tanh(kh) + kh * (1.0d0 - tanh(kh)**2)) ! derivative
      k_old = k
      k = k - f / df

      ! Check convergence
      if (abs(k - k_old) < tol * abs(k)) then
         kh = k * h
         exit
      end if
   end do

   C = omega / k
   arg = min(2.0 * kh, 50.0)
   n = 0.5 + kh / sinh(arg)
   Cg = n * C

end subroutine disper_nr
!
end module
