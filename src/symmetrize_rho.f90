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
  USE symm_base,        ONLY : nsym, s, ft
  USE fft_base,         ONLY : dffts
  USE cell_base,        ONLY : at
  !-- parameters ------------------------------------------------------
  implicit none
  real(dp), intent(inout) :: rho(dffts%nr1x,dffts%nr2x,dffts%nr3x)

  !-- local variables ----------------------------------------------------
  integer, allocatable :: symflag (:,:,:)
  integer :: ri(48), rj(48), rk(48), i, j, k, isym, ns, ftau(3)
  integer :: nr1, nr2, nr3
  real(DP) :: sum, ft_(3)

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

  nr1 = dffts%nr1
  nr2 = dffts%nr2
  nr3 = dffts%nr3
  do k = 1, nr3
     do j = 1, nr2
        do i = 1, nr1
           if (symflag(i,j,k) == 0) then
              sum = 0.d0
              do isym = 1, nsym
                 ! recover fractional translations
                 ft_(1) = ft(1,isym)*nr1
                 ft_(2) = ft(2,isym)*nr2
                 ft_(3) = ft(3,isym)*nr3
                 ftau(:) = nint(ft_(:))
                 call ruotaijk(s(1,1,isym), ftau, i, j, k, nr1, nr2, nr3, ri(isym), rj(isym), rk(isym))
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
  USE scatter_mod,      ONLY : gather_grid, scatter_grid
  !-- parameters ------------------------------------------------------
  implicit none
  real(dp), intent(inout) :: rho(dffts%nnr)

  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: aux(:)

  ! if no symmetries, return
  if (nsym <= 1) return

  allocate( aux(dffts%nr1x*dffts%nr2x*dffts%nr3x) )
  call gather_grid(dffts, rho, aux)
  if ( me_pool == 0 ) call symmetrize_rho_s(aux)
  call scatter_grid(dffts, aux, rho)

  deallocate(aux)
  return
END SUBROUTINE psymmetrize_rho_s
