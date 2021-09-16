!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE qesgipaw_reset_module
  !
  ! Auto-generated code: don't edit or at least don't commit changes
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes_1.0
  !
  USE qesgipaw_types_module
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC qesgipaw_reset
  !
  INTERFACE qesgipaw_reset
    MODULE PROCEDURE qes_reset_gipaw
    MODULE PROCEDURE qes_reset_general_info
    MODULE PROCEDURE qes_reset_format
    MODULE PROCEDURE qes_reset_creator
    MODULE PROCEDURE qes_reset_created
    MODULE PROCEDURE qes_reset_parallel_info
    MODULE PROCEDURE qes_reset_closed
    MODULE PROCEDURE qes_reset_inputGIPAW
    MODULE PROCEDURE qes_reset_job
    MODULE PROCEDURE qes_reset_files
    MODULE PROCEDURE qes_reset_scf_gipaw
    MODULE PROCEDURE qes_reset_nmr
    MODULE PROCEDURE qes_reset_efg
    MODULE PROCEDURE qes_reset_hfi
  END INTERFACE qesgipaw_reset
  !
  CONTAINS
  !
  !
  SUBROUTINE qes_reset_gipaw(obj)
    !
    IMPLICIT NONE
    TYPE(gipaw_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%general_info_ispresent) &
      CALL qes_reset_general_info(obj%general_info)
    obj%general_info_ispresent = .FALSE.
    IF (obj%parallel_info_ispresent) &
      CALL qes_reset_parallel_info(obj%parallel_info)
    obj%parallel_info_ispresent = .FALSE.
    CALL qes_reset_inputGIPAW(obj%inputGIPAW)
    obj%status_ispresent = .FALSE.
    obj%cputime_ispresent = .FALSE.
    IF (obj%closed_ispresent) &
      CALL qes_reset_closed(obj%closed)
    obj%closed_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_gipaw
  !
  !
  SUBROUTINE qes_reset_general_info(obj)
    !
    IMPLICIT NONE
    TYPE(general_info_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_format(obj%format)
    CALL qes_reset_creator(obj%creator)
    CALL qes_reset_created(obj%created)
    !
  END SUBROUTINE qes_reset_general_info
  !
  !
  SUBROUTINE qes_reset_format(obj)
    !
    IMPLICIT NONE
    TYPE(format_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_format
  !
  !
  SUBROUTINE qes_reset_creator(obj)
    !
    IMPLICIT NONE
    TYPE(creator_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_creator
  !
  !
  SUBROUTINE qes_reset_created(obj)
    !
    IMPLICIT NONE
    TYPE(created_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_created
  !
  !
  SUBROUTINE qes_reset_parallel_info(obj)
    !
    IMPLICIT NONE
    TYPE(parallel_info_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%nimages_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_parallel_info
  !
  !
  SUBROUTINE qes_reset_closed(obj)
    !
    IMPLICIT NONE
    TYPE(closed_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_closed
  !
  !
  SUBROUTINE qes_reset_inputGIPAW(obj)
    !
    IMPLICIT NONE
    TYPE(inputGIPAW_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%control_job_ispresent) &
      CALL qes_reset_job(obj%control_job)
    obj%control_job_ispresent = .FALSE.
    IF (obj%scf_gipaw_ispresent) &
      CALL qes_reset_scf_gipaw(obj%scf_gipaw)
    obj%scf_gipaw_ispresent = .FALSE.
    IF (obj%files_ispresent) &
      CALL qes_reset_files(obj%files)
    obj%files_ispresent = .FALSE.
    IF (obj%nmr_input_ispresent) &
      CALL qes_reset_nmr(obj%nmr_input)
    obj%nmr_input_ispresent = .FALSE.
    IF (obj%efg_input_ispresent) &
      CALL qes_reset_efg(obj%efg_input)
    obj%efg_input_ispresent = .FALSE.
    IF (obj%hfi_input_ispresent) &
      CALL qes_reset_hfi(obj%hfi_input)
    obj%hfi_input_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_inputGIPAW
  !
  !
  SUBROUTINE qes_reset_job(obj)
    !
    IMPLICIT NONE
    TYPE(job_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%job_ispresent = .FALSE.
    obj%restart_mode_ispresent = .FALSE.
    obj%verbosity_ispresent = .FALSE.
    obj%max_seconds_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_job
  !
  !
  SUBROUTINE qes_reset_files(obj)
    !
    IMPLICIT NONE
    TYPE(files_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%prefix_ispresent = .FALSE.
    obj%tmp_dir_ispresent = .FALSE.
    obj%filcurr_ispresent = .FALSE.
    obj%fildfield_ispresent = .FALSE.
    obj%filnics_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_files
  !
  !
  SUBROUTINE qes_reset_scf_gipaw(obj)
    !
    IMPLICIT NONE
    TYPE(scf_gipaw_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%diagonalization_ispresent = .FALSE.
    obj%conv_threshold_ispresent = .FALSE.
    obj%q_gipaw_ispresent = .FALSE.
    obj%r_rand_ispresent = .FALSE.
    obj%spline_ps_ispresent = .FALSE.
    obj%paw_proj_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_scf_gipaw
  !
  !
  SUBROUTINE qes_reset_nmr(obj)
    !
    IMPLICIT NONE
    TYPE(nmr_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%use_nmr_macroscopic_shape_ispresent = .FALSE.
    obj%nmr_macroscopic_shape_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_nmr
  !
  !
  SUBROUTINE qes_reset_efg(obj)
    !
    IMPLICIT NONE
    TYPE(efg_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%efg_q_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_efg
  !
  !
  SUBROUTINE qes_reset_hfi(obj)
    !
    IMPLICIT NONE
    TYPE(hfi_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%hfi_output_unit_ispresent = .FALSE.
    obj%hfi_nuclear_g_tensor_ispresent = .FALSE.
    obj%core_relax_method_ispresent = .FALSE.
    obj%hfi_via_reconstruction_only_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_hfi
  !
  !
END MODULE qesgipaw_reset_module