!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define __ARY

!-----------------------------------------------------------------------
SUBROUTINE knight_shift
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the knight_shift according to PRB 76, 165122 (2007)
  ! 
  USE kinds,                  ONLY : dp
#ifdef __ARY
  USE fft_base,               ONLY : dfftp
  USE scf,                    ONLY : rho, rho_core, rhog_core
  USE io_global,              ONLY : stdout
#endif

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE
#ifdef __ARY
  !integer :: s_maj, s_min
  ! the spin density
  real(dp), allocatable :: spin_den(:)
  real(dp), allocatable :: vxc(:,:)
  real(DP) :: etxc, vtxc

  integer :: ir
  real(DP) :: sumrho, sumvxc

  call start_clock('knight_shift')

  ! select majority and minority spin components
  !call select_spin(s_min, s_maj)

  !--------------------------------------------------------------------
  ! 
  !--------------------------------------------------------------------

  write(stdout,'(5X,A)') 'Computing the Knight shifts'
  write(stdout,*)

  allocate(vxc(dfftp%nnr, 2), spin_den(dfftp%nnr))

  ! recomputing the XC potentials, without Hartree contribution
  vxc(:,:) = 0.D0
  call v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )

  ! 
  spin_den(:) = abs(rho%of_r(:,1) - rho%of_r(:,2))

  sumrho = 0.D0
  sumvxc = 0.D0
  do ir = 1, dfftp%nnr

    sumrho = sumrho + spin_den(ir)
    sumvxc = sumvxc + abs(vxc(ir,1) - vxc(ir,2))

  end do

  write(stdout,'(5X,F12.6,F12.6)') sumrho / dfftp%nnr, sumvxc / dfftp%nnr
  write(stdout,*)

  ! releasing used memory
  deallocate(vxc)

  call stop_clock('knight_shift')

#endif

END SUBROUTINE knight_shift

