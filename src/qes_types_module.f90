!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE qesgipaw_types_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes_1.0 
  !
  USE kinds, only: DP
  !
  IMPLICIT NONE
  !
  TYPE :: format_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: NAME
    CHARACTER(len=256) :: VERSION
    !
    CHARACTER(len=256) :: format
    !
  END TYPE format_type
  !
  TYPE :: creator_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: NAME
    CHARACTER(len=256) :: VERSION
    !
    CHARACTER(len=256) :: creator
    !
  END TYPE creator_type
  !
  TYPE :: created_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: DATE
    CHARACTER(len=256) :: TIME
    !
    CHARACTER(len=256) :: created
    !
  END TYPE created_type
  !
  TYPE :: closed_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: DATE
    CHARACTER(len=256) :: TIME
    !
    CHARACTER(len=256) :: closed
    !
  END TYPE closed_type
  !
  TYPE :: general_info_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(format_type) :: format
    TYPE(creator_type) :: creator
    TYPE(created_type) :: created
    CHARACTER(len=256) :: job
    !
  END TYPE general_info_type
  !
  TYPE :: parallel_info_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nprocs
    INTEGER :: nthreads
    INTEGER :: ntasks
    INTEGER :: nbgrp
    INTEGER :: npool
    INTEGER :: ndiag
    LOGICAL  :: nimages_ispresent = .FALSE.
    INTEGER :: nimages
    !
  END TYPE parallel_info_type
  !
  TYPE :: job_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: job_ispresent = .FALSE.
    CHARACTER(len=256) :: job
    LOGICAL  :: restart_mode_ispresent = .FALSE.
    CHARACTER(len=256) :: restart_mode
    LOGICAL  :: verbosity_ispresent = .FALSE.
    CHARACTER(len=256) :: verbosity
    LOGICAL  :: max_seconds_ispresent = .FALSE.
    REAL(DP) :: max_seconds
    !
  END TYPE job_type
  !
  TYPE :: files_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: prefix_ispresent = .FALSE.
    CHARACTER(len=256) :: prefix
    LOGICAL  :: tmp_dir_ispresent = .FALSE.
    CHARACTER(len=256) :: tmp_dir
    LOGICAL  :: filcurr_ispresent = .FALSE.
    CHARACTER(len=256) :: filcurr
    LOGICAL  :: fildfield_ispresent = .FALSE.
    CHARACTER(len=256) :: fildfield
    LOGICAL  :: filnics_ispresent = .FALSE.
    CHARACTER(len=256) :: filnics
    !
  END TYPE files_type
  !
  TYPE :: scf_gipaw_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: diagonalization_ispresent = .FALSE.
    CHARACTER(len=256) :: diagonalization
    LOGICAL  :: conv_threshold_ispresent = .FALSE.
    REAL(DP) :: conv_threshold
    LOGICAL  :: q_gipaw_ispresent = .FALSE.
    REAL(DP) :: q_gipaw
    LOGICAL  :: r_rand_ispresent = .FALSE.
    REAL(DP) :: r_rand
    LOGICAL  :: spline_ps_ispresent = .FALSE.
    LOGICAL :: spline_ps
    LOGICAL  :: paw_proj_ispresent = .FALSE.
    LOGICAL, DIMENSION(:), ALLOCATABLE :: paw_proj
    !
  END TYPE scf_gipaw_type
  !
  TYPE :: nmr_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: use_nmr_macroscopic_shape_ispresent = .FALSE.
    LOGICAL :: use_nmr_macroscopic_shape
    LOGICAL  :: nmr_macroscopic_shape_ispresent = .FALSE.
    REAL(DP), DIMENSION(9) :: nmr_macroscopic_shape
    !
  END TYPE nmr_type
  !
  TYPE :: efg_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: efg_q_ispresent = .FALSE.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: efg_q
    !
  END TYPE efg_type
  !
  TYPE :: hfi_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: hfi_output_unit_ispresent = .FALSE.
    CHARACTER(len=256) :: hfi_output_unit
    LOGICAL  :: hfi_nuclear_g_tensor_ispresent = .FALSE.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: hfi_nuclear_g_tensor
    LOGICAL  :: core_relax_method_ispresent = .FALSE.
    INTEGER :: core_relax_method
    LOGICAL  :: hfi_via_reconstruction_only_ispresent = .FALSE.
    LOGICAL :: hfi_via_reconstruction_only
    !
  END TYPE hfi_type
  !
  TYPE :: inputGIPAW_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: control_job_ispresent = .FALSE.
    TYPE(job_type) :: control_job
    LOGICAL  :: scf_gipaw_ispresent = .FALSE.
    TYPE(scf_gipaw_type) :: scf_gipaw
    LOGICAL  :: files_ispresent = .FALSE.
    TYPE(files_type) :: files
    LOGICAL  :: nmr_input_ispresent = .FALSE.
    TYPE(nmr_type) :: nmr_input
    LOGICAL  :: efg_input_ispresent = .FALSE.
    TYPE(efg_type) :: efg_input
    LOGICAL  :: hfi_input_ispresent = .FALSE.
    TYPE(hfi_type) :: hfi_input
    !
  END TYPE inputGIPAW_type
  !
  TYPE :: gipaw_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: general_info_ispresent = .FALSE.
    TYPE(general_info_type) :: general_info
    LOGICAL  :: parallel_info_ispresent = .FALSE.
    TYPE(parallel_info_type) :: parallel_info
    TYPE(inputGIPAW_type) :: inputGIPAW
    LOGICAL  :: status_ispresent = .FALSE.
    INTEGER :: status
    LOGICAL  :: cputime_ispresent = .FALSE.
    INTEGER :: cputime
    LOGICAL  :: closed_ispresent = .FALSE.
    TYPE(closed_type) :: closed
    !
  END TYPE gipaw_type
  !
  !
END MODULE qesgipaw_types_module