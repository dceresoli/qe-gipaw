!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qesgipaw_bcast_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes_1.0
  !
  USE qesgipaw_types_module
  USE io_global, ONLY : ionode
  USE mp, ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  PUBLIC qesgipaw_bcast
  !
  INTERFACE qesgipaw_bcast
    MODULE PROCEDURE qes_bcast_gipaw
    MODULE PROCEDURE qes_bcast_general_info
    MODULE PROCEDURE qes_bcast_format
    MODULE PROCEDURE qes_bcast_creator
    MODULE PROCEDURE qes_bcast_created
    MODULE PROCEDURE qes_bcast_parallel_info
    MODULE PROCEDURE qes_bcast_closed
    MODULE PROCEDURE qes_bcast_inputGIPAW
    MODULE PROCEDURE qes_bcast_job
    MODULE PROCEDURE qes_bcast_files
    MODULE PROCEDURE qes_bcast_scf_gipaw
    MODULE PROCEDURE qes_bcast_nmr
    MODULE PROCEDURE qes_bcast_efg
    MODULE PROCEDURE qes_bcast_hfi
  END INTERFACE qesgipaw_bcast
  !
  CONTAINS
  !
  !
  SUBROUTINE qes_bcast_gipaw(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(gipaw_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%general_info_ispresent, ionode_id, comm)
    IF (obj%general_info_ispresent) &
      CALL qes_bcast_general_info(obj%general_info, ionode_id, comm)
    CALL mp_bcast(obj%parallel_info_ispresent, ionode_id, comm)
    IF (obj%parallel_info_ispresent) &
      CALL qes_bcast_parallel_info(obj%parallel_info, ionode_id, comm)
    CALL qes_bcast_inputGIPAW(obj%inputGIPAW, ionode_id, comm)
    CALL mp_bcast(obj%status_ispresent, ionode_id, comm)
    IF (obj%status_ispresent) &
      CALL mp_bcast(obj%status, ionode_id, comm)
    CALL mp_bcast(obj%cputime_ispresent, ionode_id, comm)
    IF (obj%cputime_ispresent) &
      CALL mp_bcast(obj%cputime, ionode_id, comm)
    CALL mp_bcast(obj%closed_ispresent, ionode_id, comm)
    IF (obj%closed_ispresent) &
      CALL qes_bcast_closed(obj%closed, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_gipaw
  !
  !
  SUBROUTINE qes_bcast_general_info(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(general_info_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_format(obj%format, ionode_id, comm)
    CALL qes_bcast_creator(obj%creator, ionode_id, comm)
    CALL qes_bcast_created(obj%created, ionode_id, comm)
    CALL mp_bcast(obj%job, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_general_info
  !
  !
  SUBROUTINE qes_bcast_format(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(format_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%NAME, ionode_id, comm)
    CALL mp_bcast(obj%VERSION, ionode_id, comm)
    CALL mp_bcast(obj%format, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_format
  !
  !
  SUBROUTINE qes_bcast_creator(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(creator_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%NAME, ionode_id, comm)
    CALL mp_bcast(obj%VERSION, ionode_id, comm)
    CALL mp_bcast(obj%creator, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_creator
  !
  !
  SUBROUTINE qes_bcast_created(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(created_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%DATE, ionode_id, comm)
    CALL mp_bcast(obj%TIME, ionode_id, comm)
    CALL mp_bcast(obj%created, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_created
  !
  !
  SUBROUTINE qes_bcast_parallel_info(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(parallel_info_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nprocs, ionode_id, comm)
    CALL mp_bcast(obj%nthreads, ionode_id, comm)
    CALL mp_bcast(obj%ntasks, ionode_id, comm)
    CALL mp_bcast(obj%nbgrp, ionode_id, comm)
    CALL mp_bcast(obj%npool, ionode_id, comm)
    CALL mp_bcast(obj%ndiag, ionode_id, comm)
    CALL mp_bcast(obj%nimages_ispresent, ionode_id, comm)
    IF (obj%nimages_ispresent) &
      CALL mp_bcast(obj%nimages, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_parallel_info
  !
  !
  SUBROUTINE qes_bcast_closed(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(closed_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%DATE, ionode_id, comm)
    CALL mp_bcast(obj%TIME, ionode_id, comm)
    CALL mp_bcast(obj%closed, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_closed
  !
  !
  SUBROUTINE qes_bcast_inputGIPAW(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(inputGIPAW_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%control_job_ispresent, ionode_id, comm)
    IF (obj%control_job_ispresent) &
      CALL qes_bcast_job(obj%control_job, ionode_id, comm)
    CALL mp_bcast(obj%scf_gipaw_ispresent, ionode_id, comm)
    IF (obj%scf_gipaw_ispresent) &
      CALL qes_bcast_scf_gipaw(obj%scf_gipaw, ionode_id, comm)
    CALL mp_bcast(obj%files_ispresent, ionode_id, comm)
    IF (obj%files_ispresent) &
      CALL qes_bcast_files(obj%files, ionode_id, comm)
    CALL mp_bcast(obj%nmr_input_ispresent, ionode_id, comm)
    IF (obj%nmr_input_ispresent) &
      CALL qes_bcast_nmr(obj%nmr_input, ionode_id, comm)
    CALL mp_bcast(obj%efg_input_ispresent, ionode_id, comm)
    IF (obj%efg_input_ispresent) &
      CALL qes_bcast_efg(obj%efg_input, ionode_id, comm)
    CALL mp_bcast(obj%hfi_input_ispresent, ionode_id, comm)
    IF (obj%hfi_input_ispresent) &
      CALL qes_bcast_hfi(obj%hfi_input, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_inputGIPAW
  !
  !
  SUBROUTINE qes_bcast_job(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(job_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%job_ispresent, ionode_id, comm)
    IF (obj%job_ispresent) &
      CALL mp_bcast(obj%job, ionode_id, comm)
    CALL mp_bcast(obj%restart_mode_ispresent, ionode_id, comm)
    IF (obj%restart_mode_ispresent) &
      CALL mp_bcast(obj%restart_mode, ionode_id, comm)
    CALL mp_bcast(obj%verbosity_ispresent, ionode_id, comm)
    IF (obj%verbosity_ispresent) &
      CALL mp_bcast(obj%verbosity, ionode_id, comm)
    CALL mp_bcast(obj%max_seconds_ispresent, ionode_id, comm)
    IF (obj%max_seconds_ispresent) &
      CALL mp_bcast(obj%max_seconds, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_job
  !
  !
  SUBROUTINE qes_bcast_files(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(files_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%prefix_ispresent, ionode_id, comm)
    IF (obj%prefix_ispresent) &
      CALL mp_bcast(obj%prefix, ionode_id, comm)
    CALL mp_bcast(obj%tmp_dir_ispresent, ionode_id, comm)
    IF (obj%tmp_dir_ispresent) &
      CALL mp_bcast(obj%tmp_dir, ionode_id, comm)
    CALL mp_bcast(obj%filcurr_ispresent, ionode_id, comm)
    IF (obj%filcurr_ispresent) &
      CALL mp_bcast(obj%filcurr, ionode_id, comm)
    CALL mp_bcast(obj%fildfield_ispresent, ionode_id, comm)
    IF (obj%fildfield_ispresent) &
      CALL mp_bcast(obj%fildfield, ionode_id, comm)
    CALL mp_bcast(obj%filnics_ispresent, ionode_id, comm)
    IF (obj%filnics_ispresent) &
      CALL mp_bcast(obj%filnics, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_files
  !
  !
  SUBROUTINE qes_bcast_scf_gipaw(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(scf_gipaw_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%diagonalization_ispresent, ionode_id, comm)
    IF (obj%diagonalization_ispresent) &
      CALL mp_bcast(obj%diagonalization, ionode_id, comm)
    CALL mp_bcast(obj%conv_threshold_ispresent, ionode_id, comm)
    IF (obj%conv_threshold_ispresent) &
      CALL mp_bcast(obj%conv_threshold, ionode_id, comm)
    CALL mp_bcast(obj%q_gipaw_ispresent, ionode_id, comm)
    IF (obj%q_gipaw_ispresent) &
      CALL mp_bcast(obj%q_gipaw, ionode_id, comm)
    CALL mp_bcast(obj%r_rand_ispresent, ionode_id, comm)
    IF (obj%r_rand_ispresent) &
      CALL mp_bcast(obj%r_rand, ionode_id, comm)
    CALL mp_bcast(obj%spline_ps_ispresent, ionode_id, comm)
    IF (obj%spline_ps_ispresent) &
      CALL mp_bcast(obj%spline_ps, ionode_id, comm)
    CALL mp_bcast(obj%paw_proj_ispresent, ionode_id, comm)
    IF (obj%paw_proj_ispresent) &
      CALL mp_bcast(obj%paw_proj, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_scf_gipaw
  !
  !
  SUBROUTINE qes_bcast_nmr(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(nmr_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%use_nmr_macroscopic_shape_ispresent, ionode_id, comm)
    IF (obj%use_nmr_macroscopic_shape_ispresent) &
      CALL mp_bcast(obj%use_nmr_macroscopic_shape, ionode_id, comm)
    CALL mp_bcast(obj%nmr_macroscopic_shape_ispresent, ionode_id, comm)
    IF (obj%nmr_macroscopic_shape_ispresent) &
      CALL mp_bcast(obj%nmr_macroscopic_shape, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_nmr
  !
  !
  SUBROUTINE qes_bcast_efg(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(efg_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%efg_q_ispresent, ionode_id, comm)
    IF (obj%efg_q_ispresent) &
      CALL mp_bcast(obj%efg_q, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_efg
  !
  !
  SUBROUTINE qes_bcast_hfi(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(hfi_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%hfi_output_unit_ispresent, ionode_id, comm)
    IF (obj%hfi_output_unit_ispresent) &
      CALL mp_bcast(obj%hfi_output_unit, ionode_id, comm)
    CALL mp_bcast(obj%hfi_nuclear_g_tensor_ispresent, ionode_id, comm)
    IF (obj%hfi_nuclear_g_tensor_ispresent) &
      CALL mp_bcast(obj%hfi_nuclear_g_tensor, ionode_id, comm)
    CALL mp_bcast(obj%core_relax_method_ispresent, ionode_id, comm)
    IF (obj%core_relax_method_ispresent) &
      CALL mp_bcast(obj%core_relax_method, ionode_id, comm)
    CALL mp_bcast(obj%hfi_via_reconstruction_only_ispresent, ionode_id, comm)
    IF (obj%hfi_via_reconstruction_only_ispresent) &
      CALL mp_bcast(obj%hfi_via_reconstruction_only, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_hfi
  !
  !
END MODULE qesgipaw_bcast_module