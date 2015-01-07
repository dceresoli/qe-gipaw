!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! Diagonalize and return principal axis of a 3x3 tensor
SUBROUTINE principal_axis(tens, eigs, eigv)
  USE kinds, only: dp
  IMPLICIT NONE
  real(dp), intent(in) :: tens(3,3)
  real(dp), intent(out) :: eigs(3), eigv(3,3)
  complex(dp) :: w(3,3), ev(3,3)
  real(dp) :: ei(3)
  integer :: ind(3), i, j 

  do i = 1, 3
    do j = 1, 3
       w(i,j) = cmplx(tens(i,j), 0.d0, kind=dp)
    enddo
  enddo

  w = 0.5d0 * (w + transpose(w))  
  call cdiagh(3, w, 3, ei, ev)

  ind(:) = 0
  call hpsort(3, abs(ei), ind)
  
  do i = 1, 3
     eigs(i) = ei(ind(i))
     eigv(:,i) = real(ev(:,ind(i)), kind=dp)
  enddo

  return
END SUBROUTINE principal_axis


! Diagonalize and return principal axis of a 3x3 tensor (Simpson convention)
SUBROUTINE principal_axis_simpson(tens, eigs, eigv)
  USE kinds, only: dp
  IMPLICIT NONE
  real(dp), intent(in) :: tens(3,3)
  real(dp), intent(out) :: eigs(3), eigv(3,3)
  complex(dp) :: w(3,3), ev(3,3)
  real(dp) :: ei(3), tr
  integer :: ind(3), i, j 

  do i = 1, 3
    do j = 1, 3
       w(i,j) = cmplx(tens(i,j), 0.d0, kind=dp)
    enddo
  enddo

  w = 0.5d0 * (w + transpose(w))  
  call cdiagh(3, w, 3, ei, ev)

  tr = (ei(1)+ei(2)+ei(3))/3.d0
  ind(:) = 0
  call hpsort(3, abs(ei(:)-tr), ind)
  
  do i = 1, 3
     eigs(i) = ei(ind(i))
     eigv(:,i) = real(ev(:,ind(i)), kind=dp)
  enddo

  return
END SUBROUTINE principal_axis_simpson



! Select majority and minority spin
SUBROUTINE select_spin(s_min, s_maj)
  USE kinds,        ONLY : dp
  USE scf,          ONLY : rho
  USE lsda_mod,     ONLY : nspin
#ifdef __BANDS
  USE mp_bands,     ONLY : intra_bgrp_comm
#endif
  USE mp_pools,     ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  IMPLICIT NONE
  integer, intent(out) :: s_min, s_maj
  real(dp) :: rho_diff

  rho_diff = sum(rho%of_r(:,1) - rho%of_r(:,nspin))
#ifdef __BANDS
  ! siamo sicuri di questo comunicatore?
  call mp_sum(rho_diff, intra_bgrp_comm)
#else
  call mp_sum(rho_diff, intra_pool_comm)
#endif
  if ( nspin > 1 .and. abs(rho_diff) < 1.0d-3 ) then
     call errore("select_spin", "warning, rho_diff is small", -1)
     s_maj = 1
     s_min = nspin
  else if ( rho_diff >=  0.0d0 ) then
     s_maj = 1
     s_min = nspin
  else if ( rho_diff < 0.0d0 ) then
     s_maj = nspin
     s_min = 1
  endif

END SUBROUTINE select_spin



! Trace
SUBROUTINE trace(n, a, tr)
  USE kinds,        ONLY : dp
  IMPLICIT NONE
  integer, intent(in) :: n
  real(dp), intent(in) :: a(n,n)
  real(dp), intent(out) :: tr
  integer :: i

  tr = 0.d0
  do i = 1, n
    tr = tr + a(i,i)
  enddo
END SUBROUTINE trace



! More operations with splines
SUBROUTINE radial_kinetic_energy (np, l, rdata, ydata, kin_ydata)
  USE kinds
  USE splinelib
  implicit none
  integer, intent(in) :: np, l
  real(dp), intent(in) :: rdata(np), ydata(np)
  real(dp), intent(out) :: kin_ydata(np)
  real(dp) :: d1
      
  d1 = ( ydata(2) - ydata(1) ) / ( rdata(2) - rdata(1) )
  call spline ( rdata, ydata, 0.0_dp, d1, kin_ydata )
  kin_ydata = -kin_ydata + l*(l+1) * ydata / rdata ** 2
END SUBROUTINE radial_kinetic_energy


SUBROUTINE radial_derivative (np, rdata, ydata, dydata_dr)
  USE kinds
  USE splinelib
  implicit none
  integer, intent(in) :: np
  real(dp), intent(in) :: rdata(np), ydata (np)   ! ydata passed as y * r
  real(dp), intent(out) :: dydata_dr(np)
  integer :: j
  real(dp) :: d1, tab_d2y(np)
      
  d1 = ( ydata(2) - ydata(1) ) / ( rdata(2) - rdata(1) )
  call spline ( rdata, ydata, 0.0_dp, d1, tab_d2y )
  do j = 1, np
    dydata_dr(j) = ( splint_deriv(rdata, ydata, tab_d2y, rdata(j)) - ydata(j)/rdata(j) ) / rdata(j)
  end do
END SUBROUTINE radial_derivative



!-----------------------------------------------------------------------
SUBROUTINE spherical_average(msh, r, r0, r_max, rho_g, sph)
  !-----------------------------------------------------------------------
  !
  ! ... Compute the spherical average of a density distribution
  ! ... around a point r0, up to a distance r_max
  !
  USE kinds,                 ONLY : dp
  USE constants,             ONLY : pi, tpi, fpi
  USE cell_base,             ONLY : tpiba
  USE gvect,                 ONLY : ngm, g

  !-- parameters --------------------------------------------------------
  implicit none
  integer, intent(in) :: msh
  real(dp), intent(in) :: r(msh), r0(3), r_max
  complex(dp), intent(in) :: rho_g(ngm)
  real(dp), intent(out) :: sph(msh)

  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: jl(:)
  complex(dp) :: rho0g
  real(dp) :: gg, arg
  integer :: j, ig, jmax

  sph(1:msh) = 0.d0

  ! find index corresponding to r_max
  jmax = msh
  do j = 1, msh
      if (r(j) < r_max) jmax = j
  enddo

  ! allocate Bessel function
  allocate(jl(jmax))

  ! loop over g-vectors
  do ig = 1, ngm
    gg = tpiba * sqrt(g(1,ig)*g(1,ig) + g(2,ig)*g(2,ig) + g(3,ig)*g(3,ig))

    arg = tpi * (r0(1)*g(1,ig) + r0(2)*g(2,ig) + r0(3)*g(3,ig))
    rho0g = rho_g(ig) * cmplx(cos(arg),sin(arg),kind=dp)

    call sph_bes(jmax, r, gg, 0, jl)
    sph(1:jmax) = sph(1:jmax) + real(rho0g,kind=dp) * jl(1:jmax)

  enddo

  deallocate(jl)

END SUBROUTINE spherical_average

