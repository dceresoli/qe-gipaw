!
! Copyright (C) 2008-2011 Quantum ESPRESSO and GIPAW group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE symmetrize_rho_s (rho)
  !-----------------------------------------------------------------------
  !
  !     symmetrize the smooth charge density in real space
  !
  USE kinds,            ONLY : dp
  USE symm_base,        ONLY : nsym, s, ftau
  USE fft_base,         ONLY : dffts
  !-- parameters ------------------------------------------------------
  implicit none
  real(dp), intent(inout) :: rho(dffts%nr1x,dffts%nr2x,dffts%nr3x)

  !-- local variables ----------------------------------------------------
  integer, allocatable :: symflag (:,:,:)
  integer :: ri(48), rj(48), rk(48), i, j, k, isym
  real(DP) :: sum

  ! if no symmetries, return
  if (nsym <= 1) return

  allocate (symflag(dffts%nr1x, dffts%nr2x, dffts%nr3x))    
  do k = 1, dffts%nr3
     do j = 1, dffts%nr2
        do i = 1, dffts%nr1
           symflag(i,j,k) = 0
        enddo
     enddo
  enddo

  do k = 1, dffts%nr3
     do j = 1, dffts%nr2
        do i = 1, dffts%nr1
           if (symflag(i,j,k) == 0) then
              sum = 0.d0
              do isym = 1, nsym
                 call ruotaijk(s(1,1,isym), ftau(1,isym), i, j, k, &
                      dffts%nr1, dffts%nr2, dffts%nr3, ri(isym), rj(isym), rk(isym))
                 sum = sum + rho(ri(isym),rj(isym),rk(isym))
              enddo
              sum = sum / real(nsym,kind=dp)
              ! sum contains the symmetrized charge density at point r.
              ! now fill the star of r with this sum.
              do isym = 1, nsym
                 rho(ri(isym),rj(isym),rk(isym)) = sum
                 symflag(ri(isym),rj(isym),rk(isym)) = 1
              enddo
           endif
        enddo
     enddo
  enddo

  deallocate(symflag)
  return
end subroutine symmetrize_rho_s



!-----------------------------------------------------------------------
SUBROUTINE psymmetrize_rho_s (rho)
  !-----------------------------------------------------------------------
  !
  !     symmetrize the smooth charge density in real space (parallel version)
  !
  USE kinds,            ONLY : dp  
  USE symm_base,        ONLY : nsym
  USE mp_global,        ONLY : me_pool
  USE fft_base,         ONLY : dffts
  USE fft_base,         ONLY : gather_smooth, scatter_smooth
  !-- parameters ------------------------------------------------------
  implicit none
  real(dp), intent(inout) :: rho(dffts%nnr)

  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: aux(:)

  ! if no symmetries, return
  if (nsym <= 1) return

  allocate( aux(dffts%nr1x*dffts%nr2x*dffts%nr3x) )
  call gather_smooth(rho, aux)
  if ( me_pool == 0 ) call symmetrize_rho_s(aux)
  call scatter_smooth(aux, rho)

  deallocate(aux)
  return
END SUBROUTINE psymmetrize_rho_s
