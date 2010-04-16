!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

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

