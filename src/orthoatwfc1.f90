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
  ! Is assumes that you have already called gk_sorta and init_us_2.
  !
  USE kinds,            ONLY : DP
  USE basis,            ONLY : natomwfc, swfcatom
  USE ldaU,             ONLY : Hubbard_projectors, copy_U_wfc
  USE wvfct,            ONLY : npwx
  USE klist,            ONLY : ngk
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, &
                               bec_type, becp, calbec

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik
  !-- local variables ----------------------------------------------------
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)
  INTEGER :: npw
  LOGICAL :: normalize_only

  if (Hubbard_projectors == "file") &
    call errore('orthoatwfc1', 'Hubbard_projectors == file cannot be used', 1)
 
  allocate(wfcatom(npwx,natomwfc), swfcatom(npwx,natomwfc))    
  call allocate_bec_type(nkb, natomwfc, becp)

  ! generate atomic wfcs for projection
  npw = ngk(ik)
  call atomic_wfc(ik, wfcatom)
  call calbec (npw, vkb, wfcatom, becp)
  call s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)

  if (Hubbard_projectors /= 'atomic') then
      normalize_only = (Hubbard_projectors == 'norm-atomic') 
      call ortho_swfc(npw, normalize_only, natomwfc, wfcatom, swfcatom, .FALSE.)
  endif

  CALL copy_U_wfc (swfcatom)

  call deallocate_bec_type(becp)
  deallocate(wfcatom, swfcatom)

END SUBROUTINE orthoatwfc1

