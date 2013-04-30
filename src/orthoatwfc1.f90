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
  USE basis,            ONLY : natomwfc, swfcatom
  USE klist,            ONLY : nks, xk, ngk
  USE ldaU,             ONLY : U_projection, wfcU, nwfcU, copy_U_wfc
  USE wvfct,            ONLY : npwx, npw, igk
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, &
                               bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  USE io_global,        ONLY : stdout

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  integer, intent(in) :: ik
  !-- local variables ----------------------------------------------------
  complex(dp), allocatable :: wfcatom(:,:)
  LOGICAL :: normalize_only

  allocate(wfcatom(npwx,natomwfc))    

  if (U_projection == "file") &
    call errore('orthoatwfc1', 'U_projection == file cannot be used', 1)
 
  ! generate atomic wfcs for projection
  call atomic_wfc(ik, wfcatom)
  call s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)

  if (U_projection=="atomic") goto 200

  ! orthogonalize

  normalize_only = (U_projection == "norm-atomic") 
  CALL ortho_swfc ( normalize_only, natomwfc, wfcatom, swfcatom )
  !
200 continue
  !
  CALL copy_U_wfc (swfcatom)
  !
  deallocate(wfcatom)

END SUBROUTINE orthoatwfc1
