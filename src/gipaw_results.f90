!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE gipaw_results
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the results of GIPAW calculations
  !
  USE kinds,             ONLY : dp
  
  IMPLICIT NONE
  SAVE

  ! GIPAW results, in order to dump them on the XML file
  real(dp) :: res_suscept1(3,3), res_suscept2(3,3)
  real(dp), allocatable :: res_nmr_sigma(:,:,:)
  real(dp), allocatable :: res_efg(:,:,:)
  real(dp) :: res_epr_deltag(3,3), res_epr_deltag_paratec(3,3)
  real(dp), allocatable :: res_hfi_dip(:,:,:), res_hfi_fc(:)

  PUBLIC :: allocate_gipaw_results

  CONTAINS

  SUBROUTINE allocate_gipaw_results
    USE ions_base,     ONLY : ntyp => nsp, nat
    implicit none

    allocate(res_nmr_sigma(3,3,nat))
    allocate(res_efg(3,3,nat), res_hfi_dip(3,3,nat), res_hfi_fc(nat))

    res_suscept1 = 0.d0
    res_suscept2 = 0.d0
    res_nmr_sigma = 0.d0
    res_efg = 0.d0
    res_epr_deltag = 0.d0
    res_epr_deltag_paratec = 0.d0
    res_hfi_dip = 0.d0
    res_hfi_fc = 0.d0

  END SUBROUTINE allocate_gipaw_results

END MODULE gipaw_results
