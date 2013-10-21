!
! Copyright (C) 2001-2013 Quantum-Espresso Group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! see also gipaw_routines_bands.f90 for s_spi_bands and add_vuspsi_bands
!
! TODO: Oct 21 2013 (DC): the band parallelization has to tested again
!
!-----------------------------------------------------------------------
#ifdef __BANDS
subroutine h_psiq_bands (lda, n, m, psi, hpsi, spsi, ibnd_start, ibnd_end)
#else
subroutine h_psiq (lda, n, m, psi, hpsi, spsi)
#endif
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with a m  wavefunctions  contained
  !     in psi. It first computes the bec matrix of these
  !     wavefunctions and then with the routines hus_1psi and
  !     s_psi computes for each band the required products
  !
  USE kinds,                 ONLY : dp
  USE becmod,                ONLY : becp, calbec
  USE uspp,                  ONLY : vkb, nkb
#ifdef __BANDS
  USE mp,                    ONLY : mp_sum
  USE lsda_mod,              ONLY : current_spin
  USE wvfct,                 ONLY : g2kin
  USE scf,                   ONLY : vrs
  USE noncollin_module,      ONLY : npol, noncolin
#endif
  !--------------------------------------------------------------------
  IMPLICIT NONE
  integer, intent(in) :: lda     ! the leading dimension of the array psi (npwx)
  integer, intent(in) :: n       ! the real dimension of psi (npw)
  integer, intent(in) :: m       ! the number of psi to compute (nbnd)
#ifdef __BANDS
  integer, intent(in) :: ibnd_start, ibnd_end
  integer :: ibnd
#endif
  complex(dp), intent(in) :: psi(lda,m)
  complex(dp), intent(out) :: hpsi(lda,m), spsi(lda,m)

  call start_clock ('h_psiq')

  hpsi(:,:) = (0.d0,0.d0)

#ifdef __BANDS

  ! Here we apply the kinetic energy (k+G)^2 psi
  do ibnd = ibnd_start, ibnd_end
     hpsi (1:n, ibnd) = g2kin (1:n) * psi (1:n, ibnd)
     hpsi (n+1:lda,ibnd) = (0.0_dp, 0.0_dp)
     if (noncolin) then
        hpsi(lda+1:lda+n, ibnd) = g2kin(1:n) * psi(lda+1:lda+n, ibnd)
        hpsi(lda+n+1:lda*npol, ibnd) = (0.0_dp, 0.0_dp)
     endif
  enddo
  ! Here we add the Hubbard potential times psi
  if (lda_plus_u .and. U_projection  /= "pseudo" ) then
     if (noncolin) then
        call vhpsi_nc(lda, n, ibnd_end-ibnd_start+1, psi(1,ibnd_start), hpsi(1,ibnd_start) )
     else
        call vhpsi(lda, n, ibnd_end-ibnd_start+1, psi(1,ibnd_start), hpsi(1,ibnd_start) )
     endif
  endif
  ! Here we add the local potential
  call vloc_psi_k(lda, n, ibnd_end-ibnd_start+1, psi(1,ibnd_start), &
                  vrs(1,current_spin), hpsi(1,ibnd_start) )
  ! Here we add the NL potential
  call calbec_bands(lda, n, nkb, vkb, psi, becp%k, m, ibnd_start, ibnd_end)
  call add_vuspsi_bands (lda, n, m, hpsi, ibnd_start, ibnd_end)
  ! Here we apply the overlap
  call s_psi_bands(lda, n, m, psi, spsi, ibnd_start, ibnd_end)

#else

  call h_psi(lda, n, m, psi, hpsi)
  call calbec(n, vkb, psi, becp, m)
  call s_psi(lda, n, m, psi, spsi)

#endif

  call stop_clock ('h_psiq')
  return

end subroutine

