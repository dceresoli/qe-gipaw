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
  real(dp) :: suscept1(3,3), suscept2(3,3)
  real(dp), allocatable :: nmr_sigma_tot(:,:,:)
  real(dp) :: epr_Dg_tensor(3,3), epr_Dg_tensor_paratec(3,3)
  real(dp), allocatable :: efg_efg_tot(:,:,:)
  real(dp), allocatable :: hfi_dip_tot(:,:,:), hfi_fc_tot(:)

  PUBLIC :: allocate_gipaw_results

  CONTAINS

  SUBROUTINE allocate_gipaw_results
    USE ions_base,     ONLY : ntyp => nsp, nat
    implicit none

    allocate(nmr_sigma_tot(3,3,nat))
    allocate(efg_efg_tot(3,3,nat), hfi_dip_tot(3,3,nat), hfi_fc_tot(nat))

  END SUBROUTINE allocate_gipaw_results

END MODULE gipaw_results
