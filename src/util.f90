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



! Selent majority and minority spin
SUBROUTINE select_spin(s_min, s_maj)
  USE kinds,        ONLY : dp
  USE scf,          ONLY : rho
  USE lsda_mod,     ONLY : nspin
  USE mp_global,    ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  IMPLICIT NONE
  integer, intent(out) :: s_min, s_maj
  real(dp) :: rho_diff

  rho_diff = sum(rho%of_r(:,1) - rho%of_r(:,nspin))
  call mp_sum(rho_diff, intra_pool_comm)

  if ( nspin > 1 .and. abs(rho_diff) < 1.0d-3 ) &
     call errore("select_spin", "warning, rho_diff is small", -1)

  if ( rho_diff >=  0.0d0 ) then
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
  
