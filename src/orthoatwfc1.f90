!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthoatwfc1(ik)
  !-----------------------------------------------------------------------
  !
  ! This routine is meant to orthogonalize atomic wfcs at one k-point.
  ! Is assumes that you have already called gk_sort, init_us_2 and ccalbec
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : nks, xk, ngk
  USE ldaU,             ONLY : swfcatom, U_projection, wfcU, nwfcU, copy_U_wfc
  USE wvfct,            ONLY : npwx, npw, igk
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, &
                               bec_type, becp, calbec
  USE mp_global,        ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  USE io_global,        ONLY : stdout

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  integer, intent(in) :: ik
  !-- local variables ----------------------------------------------------
  complex(dp), allocatable :: wfcatom(:,:), work(:,:), overlap(:,:)
  real(dp), allocatable :: e(:)
  complex :: temp
  integer :: i, j, k

  allocate(wfcatom(npwx,natomwfc))    

  if (U_projection == "file") &
    call errore('orthoatwfc1', 'U_projection == file cannot be used', 1)
 
  ! generate atomic wfcs for projection
  call atomic_wfc(ik, wfcatom)
  call s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)

  if (U_projection=="atomic") goto 200

  ! orthogonalize
  allocate(overlap(natomwfc,natomwfc), work(natomwfc,natomwfc))    
  allocate(e(natomwfc))    
  overlap(:,:) = (0.d0,0.d0)
  work(:,:) = (0.d0,0.d0)

  ! calculate overlap
  call zgemm('c', 'n', natomwfc, natomwfc, npw, (1.d0,0.d0), &
             wfcatom, npwx, swfcatom, npwx, (0.d0, 0.d0), overlap, natomwfc)
#ifdef __MPI
  call mp_sum(overlap, intra_pool_comm)
#endif

  if (U_projection == "norm-atomic") then
    do i = 1, natomwfc
      do j = i+1, natomwfc
        overlap(i,j) = (0.d0,0.d0)
        overlap(j,i) = (0.d0,0.d0)
      enddo
    enddo
  endif

  ! find O^-.5
  call cdiagh(natomwfc, overlap, natomwfc, e, work)
  do i = 1, natomwfc
     e(i) = 1.d0/sqrt(e(i))
  enddo

  do i = 1, natomwfc
    do j = i, natomwfc
      temp = (0.d0, 0.d0)
      do k = 1, natomwfc
         temp = temp + e(k) * work(j,k) * conjg(work(i,k))
      enddo
      overlap(i,j) = temp
      if (i /= j) overlap(j,i) = conjg(temp)
     enddo
  enddo

  ! trasform atomic orbitals O^-.5 psi
  do i = 1, npw
    work(:,1) = (0.d0,0.d0)
    call zgemv('n', natomwfc, natomwfc, (1.d0,0.d0), overlap, &
                natomwfc, swfcatom (i,1), npwx, (0.d0,0.d0), work, 1)
    call zcopy(natomwfc, work, 1, swfcatom (i,1), npwx)
  enddo
  deallocate(overlap, work, e)
        
200 continue
  !!!
  CALL copy_U_wfc ()
  !!!
  deallocate(overlap, work, e, wfcatom)

END SUBROUTINE orthoatwfc1
