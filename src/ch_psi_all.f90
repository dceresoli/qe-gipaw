!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
#ifdef __BANDS
subroutine ch_psi_all (n, h, ah, e, ik, m, ibnd_start, ibnd_end)
#else
subroutine ch_psi_all (n, h, ah, e, ik, m)
#endif
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  USE kinds,        ONLY : dp
  USE wvfct,        ONLY : npwx, nbnd
  USE uspp,         ONLY : vkb, nkb
  USE becmod,       ONLY : becp, calbec
  USE gipaw_module, ONLY : nbnd_occ, alpha_pv, evq
  USE mp_pools,     ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
#ifdef __BANDS
  USE mp_bands,     ONLY : intra_bgrp_comm
#endif
  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point
#ifdef __BANDS
  integer, intent(in) :: ibnd_start, ibnd_end
#endif
  real(DP) :: e (m)
  ! input: the eigenvalue

  complex(DP) :: h (npwx, m), ah (npwx, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  complex(DP), allocatable :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h

  call start_clock ('ch_psi')
  allocate (ps  ( nbnd , m))    
  allocate (hpsi( npwx , m))    
  allocate (spsi( npwx , m))    
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
#ifdef __BANDS
  call h_psiq_bands (npwx, n, m, h, hpsi, spsi, ibnd_start, ibnd_end)
#else
  call h_psiq (npwx, n, m, h, hpsi, spsi)
#endif
  call start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  ah = (0.d0,0.d0)
#ifdef __BANDS
  do ibnd = ibnd_start, ibnd_end
#else
  do ibnd = 1, m
#endif
     do ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     enddo
  enddo
  !
  !   Here we compute the projector in the valence band
  !
  ikq = ik
  ps (:,:) = (0.d0, 0.d0)
#ifdef __BANDS
  call zgemm ('C', 'N', nbnd_occ (ikq) , ibnd_end-ibnd_start+1, n, (1.d0, 0.d0) , evq, &
       npwx, spsi(1,ibnd_start), npwx, (0.d0, 0.d0) , ps(1,ibnd_start), nbnd)
  ps (:,ibnd_start:ibnd_end) = ps(:,ibnd_start:ibnd_end) * alpha_pv
#else
  call zgemm ('C', 'N', nbnd_occ (ikq) , m, n, (1.d0, 0.d0) , evq, &
       npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
  ps (:,:) = ps(:,:) * alpha_pv
#endif

#ifdef __MPI
#ifdef __BANDS
  call mp_sum ( ps, intra_bgrp_comm )
#else
  call mp_sum ( ps, intra_pool_comm )
#endif
#endif

  hpsi (:,:) = (0.d0, 0.d0)
#ifdef __BANDS
  call zgemm ('N', 'N', n, ibnd_end-ibnd_start+1, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
       npwx, ps(1,ibnd_start), nbnd, (1.d0, 0.d0) , hpsi(1,ibnd_start), npwx)
  spsi(:,ibnd_start:ibnd_end) = hpsi(:,ibnd_start:ibnd_end)
#else
  call zgemm ('N', 'N', n, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
       npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
  spsi(:,:) = hpsi(:,:)
#endif
  !
  !    And apply S again
  !
#ifdef __BANDS
  call calbec_bands(npwx, n, nkb, vkb, hpsi, becp%k, m, ibnd_start, ibnd_end)
  call s_psi_bands (npwx, n, m, hpsi, spsi, ibnd_start, ibnd_end)
#else
  call calbec(n, vkb, hpsi, becp, m)
  call s_psi (npwx, n, m, hpsi, spsi)
#endif

#ifdef __BANDS
  do ibnd = ibnd_start, ibnd_end 
#else
  do ibnd = 1, m
#endif
     do ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
     enddo
  enddo

  deallocate (spsi)
  deallocate (hpsi)
  deallocate (ps)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
end subroutine ch_psi_all
