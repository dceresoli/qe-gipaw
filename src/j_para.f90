! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE j_para(fact, psi_n, psi_m, ik, q, j)
  !-----------------------------------------------------------------------
  !
  ! ... Compute the paramgnetic current between two states:
  ! ... j_para(r') = (1/2) fact <psi_n| { p_k|r'><r'| + |r'><r'|p_{k+q} } |psi_m>
  ! ... the result is added to j and is returned in real space
  ! ...
  ! ... In the USPP case, current is computed on the smooth grid.
  !  
  USE kinds,                  ONLY : DP
  USE klist,                  ONLY : xk
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg
  USE smooth_grid_dimensions, ONLY : nrxxs
  USE gvecs,                  ONLY : nls
  USE gvect,                  ONLY : g
  USE cell_base,              ONLY : tpiba
  USE gipaw_module,           ONLY : nbnd_occ
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : invfft
  USE mp_image_global_module,      ONLY : inter_image_comm, nimage, my_image_id
  USE mp,                          ONLY : mp_sum 
  USE mp_global,                   ONLY : inter_bgrp_comm
  USE gipaw_module,                ONLY : ibnd_start, ibnd_end
  
  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik               ! k-point
  REAL(DP), INTENT(IN) :: fact            ! multiplication factor
  REAL(DP), INTENT(IN) :: q(3)
  COMPLEX(DP), INTENT(IN) :: psi_n(npwx,nbnd), psi_m(npwx,nbnd)
  REAL(DP), INTENT(INOUT) :: j(nrxxs,3)

  !-- local variables ----------------------------------------------------
  REAL(DP), allocatable :: j_new(:,:) 
  COMPLEX(DP), allocatable :: p_psic(:), psic(:), aux(:)
  REAL(DP) :: gk
  INTEGER :: ig, ipol, ibnd, ierr

  call start_clock('j_para')

  ! allocate real space wavefunctions
  allocate(p_psic(nrxxs), psic(nrxxs), aux(npwx), j_new(nrxxs,3))

  j_new = 0.0d0  
  ! loop over cartesian component
  do ipol = 1, 3
  
    ! loop over bands
    !do ibnd = 1, nbnd_occ(ik)
    do ibnd = ibnd_start, ibnd_end
    
      ! apply p_k on the left
      do ig = 1, npw
        gk = xk(ipol,ik) + g(ipol,igk(ig))
        aux(ig) = gk * tpiba * psi_n(ig,ibnd)
      enddo
     
      ! transform to real space
      p_psic(:) = (0.d0,0.d0)
      p_psic(nls(igk(1:npw))) = aux(1:npw)
      CALL invfft ('Wave', p_psic, dffts)

      psic(:) = (0.d0,0.d0)
      psic(nls(igk(1:npw))) = psi_m(1:npw,ibnd)
      CALL invfft ('Wave', psic, dffts)

      ! add to the current
      j_new(1:nrxxs,ipol) = j_new(1:nrxxs,ipol) + 0.5d0 * fact * wg(ibnd,ik) * &
                        aimag(conjg(p_psic(1:nrxxs)) * psic(1:nrxxs))

      ! apply p_{k+q} on the right
      do ig = 1, npw
        gk = xk(ipol,ik) + g(ipol,igk(ig)) + q(ipol)
        aux(ig) = gk * tpiba * psi_m(ig,ibnd)
      enddo
     
      ! transform to real space
      p_psic(:) = (0.d0,0.d0)
      p_psic(nls(igk(1:npw))) = aux(1:npw)
      CALL invfft ('Wave', p_psic, dffts)

      psic(:) = (0.d0,0.d0)
      psic(nls(igk(1:npw))) = psi_n(1:npw,ibnd)
      CALL invfft ('Wave', psic, dffts)

      ! add to the current
      j_new(1:nrxxs,ipol) = j_new(1:nrxxs,ipol) + 0.5d0 * fact * wg(ibnd,ik) * &
                        aimag(conjg(psic(1:nrxxs)) * p_psic(1:nrxxs))

    enddo ! ibnd

    !CALL MPI_ALLGATHERV( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, j_new(1,ipol),  rcount, displ, MPI_DOUBLE_PRECISION, inter_image_comm, ierr)
#ifdef __PARA 
#ifdef __BANDS
    call mp_sum(j_new(1:nrxxs,ipol),inter_bgrp_comm)
#endif
#endif
  enddo ! ipol
  ! free memory
  j = j+j_new
  deallocate(p_psic, psic, aux, j_new)

  call stop_clock('j_para')

END SUBROUTINE j_para

 
