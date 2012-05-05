!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! see also gipaw_routines_bands.f90 for s_spi_bands and add_vuspsi_bands
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
  USE fft_base,              ONLY : dffts
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE gvecs,                 ONLY : nls
  USE lsda_mod,              ONLY : current_spin
  USE wvfct,                 ONLY : igk, g2kin
  USE uspp,                  ONLY : vkb, nkb
  USE scf,                   ONLY : vrs
  USE wavefunctions_module,  ONLY : psic
  USE becmod,                ONLY : becp, calbec
#ifdef __BANDS
  USE mp,                    ONLY : mp_sum
#endif
  implicit none
  !
  !     Here the local variables
  !
#ifdef __BANDS
  integer :: ibnd_start, ibnd_end, ibnd
#else
  integer :: ibnd
#endif
  ! counter on bands

  integer :: lda, n, m
  ! input: the leading dimension of the array psi
  ! input: the real dimension of psi
  ! input: the number of psi to compute
  integer :: j
  ! do loop index

  complex(DP) :: psi (lda, m), hpsi (lda, m), spsi (lda, m)
  ! input: the functions where to apply H and S
  ! output: H times psi
  ! output: S times psi (Us PP's only)


  call start_clock ('h_psiq')
  call start_clock ('init')

#ifdef __BANDS
  call calbec_bands (lda, n, nkb, vkb, psi, becp%k, m, ibnd_start, ibnd_end)
#else
  call calbec (n, vkb, psi, becp, m)
#endif
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  hpsi(:,:) = (0.d0,0.d0)
#ifdef __BANDS
  do ibnd = ibnd_start, ibnd_end
#else
  do ibnd = 1, m
#endif
     do j = 1, n
        hpsi (j, ibnd) = g2kin (j) * psi (j, ibnd)
     enddo
  enddo
  call stop_clock ('init')
  !
  ! the local potential V_Loc psi. First the psi in real space
  !

#ifdef __BANDS
  do ibnd = ibnd_start, ibnd_end
#else
  do ibnd = 1, m
#endif
     call start_clock ('firstfft')
     psic(:) = (0.d0, 0.d0)
     psic (nls(igk(1:n))) = psi(1:n, ibnd)
     CALL invfft ('Wave', psic, dffts)
     call stop_clock ('firstfft')
     !
     !   and then the product with the potential vrs = (vltot+vr) on the smoo
     !
     psic (1:dffts%nnr) = psic (1:dffts%nnr) * vrs (1:dffts%nnr, current_spin)
     !
     !   back to reciprocal space
     !
     call start_clock ('secondfft')
     CALL fwfft ('Wave', psic, dffts)
     !
     !   addition to the total product
     !
     hpsi (1:n, ibnd) = hpsi (1:n, ibnd) + psic (nls(igk(1:n)))
     call stop_clock ('secondfft')
  enddo
  !
  !  Here the product with the non local potential V_NL psi
  !
#ifdef __BANDS
  call add_vuspsi_bands (lda, n, m, hpsi, ibnd_start, ibnd_end)
  call s_psi_bands (lda, n, m, psi, spsi, ibnd_start, ibnd_end)
#else
  call add_vuspsi (lda, n, m, hpsi)
  call s_psi (lda, n, m, psi, spsi)
#endif

  call stop_clock ('h_psiq')
  return

end subroutine
