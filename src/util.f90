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


! transform rho from (rho,zeta) to (up,down)
SUBROUTINE get_rho_up_down
  USE kinds,        ONLY : dp
  USE scf,          ONLY : rho
  USE lsda_mod,     ONLY : nspin
  USE fft_base,     ONLY : dfftp
  USE io_global,    ONLY : stdout
  IMPLICIT NONE
  real(dp) :: r, z
  integer :: i

  if (nspin == 1) return

  do i = 1, dfftp%nnr
     r = rho%of_r(i,1)
     z = rho%of_r(i,2)
     rho%of_r(i,1) = 0.5d0*(r+z)
     rho%of_r(i,2) = 0.5d0*(r-z)
  enddo
  write(stdout,'(/,5X,''(RHO,ZETA) => (RHO_UP,RHO_DOWN)'',/)')

END SUBROUTINE get_rho_up_down


! transform rho from (up,down) to (rho,zeta)
SUBROUTINE get_rho_zeta
  USE kinds,        ONLY : dp
  USE scf,          ONLY : rho
  USE lsda_mod,     ONLY : nspin
  USE fft_base,     ONLY : dfftp
  USE io_global,    ONLY : stdout
  IMPLICIT NONE
  real(dp) :: up, dw
  integer :: i

  if (nspin == 1) return

  do i = 1, dfftp%nnr
     up = rho%of_r(i,1)
     dw = rho%of_r(i,2)
     rho%of_r(i,1) = up + dw
     rho%of_r(i,2) = up - dw
  enddo
  write(stdout,'(/,5X,''(RHO_UP,RHO_DOWN) => (RHO,ZETA)'',/)')

END SUBROUTINE get_rho_zeta



! Select majority and minority spin
! starting from QE-6.4, rho has to be transformed back to (up,down)
! before calling this routine
SUBROUTINE select_spin(s_min, s_maj)
  USE kinds,        ONLY : dp
  USE scf,          ONLY : rho
  USE lsda_mod,     ONLY : nspin
#ifdef __BANDS
  USE mp_bands,     ONLY : intra_bgrp_comm
#endif
  USE mp_pools,     ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE cell_base,    ONLY : omega
  USE io_global,    ONLY : stdout
  USE fft_base,     ONLY : dfftp
  IMPLICIT NONE
  integer, intent(out) :: s_min, s_maj
  integer :: nrxx
  real(dp) :: rho_diff

  if (nspin == 1) then
      s_maj = 1
      s_min = 1
      return
  endif

  ! calculate spin magnetization
  nrxx = dfftp%nnr
  rho_diff = sum(rho%of_r(1:nrxx,1) - rho%of_r(1:nrxx,nspin))
#ifdef __BANDS
  ! siamo sicuri di questo comunicatore?
  call mp_sum(rho_diff, intra_bgrp_comm)
#else
  call mp_sum(rho_diff, intra_pool_comm)
#endif
  rho_diff = rho_diff * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)

  if ( abs(rho_diff) < 1.0d-3 ) then
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
  write(stdout,'(5X,''select_spin: s_maj='',I1,'' s_min='',I1,'' rho_diff='',F12.6)') s_maj, s_min, rho_diff
  write(stdout,*)

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


!----------------------------------------------------------------------
subroutine ruotaijk (ss, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk)
  !----------------------------------------------------------------------
  !
  !    This routine computes the rotated of the point i,j,k throught
  !    the symmetry (s,f). Then it computes the equivalent point
  !    on the original mesh
  !
  !
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: ss (3, 3), ftau (3), i, j, k, nr1, nr2, nr3, ri, rj, rk
  ! input: the rotation matrix
  ! input: the fractionary translation
  !   !   input: the point to rotate

  ! /
  !   !   input: the dimension of the mesh

  ! /
  !  !  output: the rotated point

  !/
  !
  ri = ss (1, 1) * (i - 1) + ss (2, 1) * (j - 1) + ss (3, 1) &
       * (k - 1) - ftau (1)
  ri = mod (ri, nr1) + 1
  if (ri.lt.1) ri = ri + nr1
  rj = ss (1, 2) * (i - 1) + ss (2, 2) * (j - 1) + ss (3, 2) &
       * (k - 1) - ftau (2)
  rj = mod (rj, nr2) + 1
  if (rj.lt.1) rj = rj + nr2
  rk = ss (1, 3) * (i - 1) + ss (2, 3) * (j - 1) + ss (3, 3) &
       * (k - 1) - ftau (3)
  rk = mod (rk, nr3) + 1
  if (rk.lt.1) rk = rk + nr3

  return
end subroutine ruotaijk

