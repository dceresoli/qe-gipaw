!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qesgipaw_write_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes_1.0
  !
  USE qesgipaw_types_module
  USE FoX_wxml
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
  INTERFACE qesgipaw_write
    MODULE PROCEDURE qes_write_gipaw
    MODULE PROCEDURE qes_write_general_info
    MODULE PROCEDURE qes_write_format
    MODULE PROCEDURE qes_write_creator
    MODULE PROCEDURE qes_write_created
    MODULE PROCEDURE qes_write_parallel_info
    MODULE PROCEDURE qes_write_closed
    MODULE PROCEDURE qes_write_inputGIPAW
    MODULE PROCEDURE qes_write_job
    MODULE PROCEDURE qes_write_files
    MODULE PROCEDURE qes_write_scf_gipaw
    MODULE PROCEDURE qes_write_nmr
    MODULE PROCEDURE qes_write_efg
    MODULE PROCEDURE qes_write_hfi
  END INTERFACE qesgipaw_write
  !
  CONTAINS
  !
  
   SUBROUTINE qes_write_gipaw(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(gipaw_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%general_info_ispresent) THEN
        CALL qes_write_general_info (xp, obj%general_info)
     END IF
     IF (obj%parallel_info_ispresent) THEN
        CALL qes_write_parallel_info (xp, obj%parallel_info)
     END IF
     CALL qes_write_inputGIPAW (xp, obj%inputGIPAW)
     IF (obj%status_ispresent) THEN
        CALL xml_NewElement(xp, "status")
           CALL xml_addCharacters(xp, obj%status)
        CALL xml_EndElement(xp, "status")
     END IF
     IF (obj%cputime_ispresent) THEN
        CALL xml_NewElement(xp, "cputime")
           CALL xml_addCharacters(xp, obj%cputime)
        CALL xml_EndElement(xp, "cputime")
     END IF
     IF (obj%closed_ispresent) THEN
        CALL qes_write_closed (xp, obj%closed)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_gipaw

   SUBROUTINE qes_write_general_info(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(general_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_format (xp, obj%format)
     CALL qes_write_creator (xp, obj%creator)
     CALL qes_write_created (xp, obj%created)
     CALL xml_NewElement(xp, 'job')
        CALL xml_addCharacters(xp, TRIM(obj%job))
     CALL xml_EndElement(xp, 'job')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_general_info

   SUBROUTINE qes_write_format(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(format_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'NAME', TRIM(obj%NAME) )
     CALL xml_addAttribute(xp, 'VERSION', TRIM(obj%VERSION) )
        CALL xml_AddCharacters(xp, TRIM(obj%format))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_format

   SUBROUTINE qes_write_creator(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(creator_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'NAME', TRIM(obj%NAME) )
     CALL xml_addAttribute(xp, 'VERSION', TRIM(obj%VERSION) )
        CALL xml_AddCharacters(xp, TRIM(obj%creator))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_creator

   SUBROUTINE qes_write_created(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(created_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'DATE', TRIM(obj%DATE) )
     CALL xml_addAttribute(xp, 'TIME', TRIM(obj%TIME) )
        CALL xml_AddCharacters(xp, TRIM(obj%created))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_created

   SUBROUTINE qes_write_parallel_info(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(parallel_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'nprocs')
        CALL xml_addCharacters(xp, obj%nprocs)
     CALL xml_EndElement(xp, 'nprocs')
     CALL xml_NewElement(xp, 'nthreads')
        CALL xml_addCharacters(xp, obj%nthreads)
     CALL xml_EndElement(xp, 'nthreads')
     CALL xml_NewElement(xp, 'ntasks')
        CALL xml_addCharacters(xp, obj%ntasks)
     CALL xml_EndElement(xp, 'ntasks')
     CALL xml_NewElement(xp, 'nbgrp')
        CALL xml_addCharacters(xp, obj%nbgrp)
     CALL xml_EndElement(xp, 'nbgrp')
     CALL xml_NewElement(xp, 'npool')
        CALL xml_addCharacters(xp, obj%npool)
     CALL xml_EndElement(xp, 'npool')
     CALL xml_NewElement(xp, 'ndiag')
        CALL xml_addCharacters(xp, obj%ndiag)
     CALL xml_EndElement(xp, 'ndiag')
     IF (obj%nimages_ispresent) THEN
        CALL xml_NewElement(xp, "nimages")
           CALL xml_addCharacters(xp, obj%nimages)
        CALL xml_EndElement(xp, "nimages")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_parallel_info

   SUBROUTINE qes_write_closed(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(closed_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'DATE', TRIM(obj%DATE) )
     CALL xml_addAttribute(xp, 'TIME', TRIM(obj%TIME) )
        CALL xml_AddCharacters(xp, TRIM(obj%closed))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_closed

   SUBROUTINE qes_write_inputGIPAW(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(inputGIPAW_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%control_job_ispresent) THEN
        CALL qes_write_job (xp, obj%control_job)
     END IF
     IF (obj%scf_gipaw_ispresent) THEN
        CALL qes_write_scf_gipaw (xp, obj%scf_gipaw)
     END IF
     IF (obj%files_ispresent) THEN
        CALL qes_write_files (xp, obj%files)
     END IF
     IF (obj%nmr_input_ispresent) THEN
        CALL qes_write_nmr (xp, obj%nmr_input)
     END IF
     IF (obj%efg_input_ispresent) THEN
        CALL qes_write_efg (xp, obj%efg_input)
     END IF
     IF (obj%hfi_input_ispresent) THEN
        CALL qes_write_hfi (xp, obj%hfi_input)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_inputGIPAW

   SUBROUTINE qes_write_job(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(job_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%job_ispresent) THEN
        CALL xml_NewElement(xp, "job")
           CALL xml_addCharacters(xp, TRIM(obj%job))
        CALL xml_EndElement(xp, "job")
     END IF
     IF (obj%restart_mode_ispresent) THEN
        CALL xml_NewElement(xp, "restart_mode")
           CALL xml_addCharacters(xp, TRIM(obj%restart_mode))
        CALL xml_EndElement(xp, "restart_mode")
     END IF
     IF (obj%verbosity_ispresent) THEN
        CALL xml_NewElement(xp, "verbosity")
           CALL xml_addCharacters(xp, TRIM(obj%verbosity))
        CALL xml_EndElement(xp, "verbosity")
     END IF
     IF (obj%max_seconds_ispresent) THEN
        CALL xml_NewElement(xp, "max_seconds")
           CALL xml_addCharacters(xp, obj%max_seconds, fmt='s16')
        CALL xml_EndElement(xp, "max_seconds")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_job

   SUBROUTINE qes_write_files(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(files_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%prefix_ispresent) THEN
        CALL xml_NewElement(xp, "prefix")
           CALL xml_addCharacters(xp, TRIM(obj%prefix))
        CALL xml_EndElement(xp, "prefix")
     END IF
     IF (obj%tmp_dir_ispresent) THEN
        CALL xml_NewElement(xp, "tmp_dir")
           CALL xml_addCharacters(xp, TRIM(obj%tmp_dir))
        CALL xml_EndElement(xp, "tmp_dir")
     END IF
     IF (obj%filcurr_ispresent) THEN
        CALL xml_NewElement(xp, "filcurr")
           CALL xml_addCharacters(xp, TRIM(obj%filcurr))
        CALL xml_EndElement(xp, "filcurr")
     END IF
     IF (obj%fildfield_ispresent) THEN
        CALL xml_NewElement(xp, "fildfield")
           CALL xml_addCharacters(xp, TRIM(obj%fildfield))
        CALL xml_EndElement(xp, "fildfield")
     END IF
     IF (obj%filnics_ispresent) THEN
        CALL xml_NewElement(xp, "filnics")
           CALL xml_addCharacters(xp, TRIM(obj%filnics))
        CALL xml_EndElement(xp, "filnics")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_files

   SUBROUTINE qes_write_scf_gipaw(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(scf_gipaw_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%diagonalization_ispresent) THEN
        CALL xml_NewElement(xp, "diagonalization")
           CALL xml_addCharacters(xp, TRIM(obj%diagonalization))
        CALL xml_EndElement(xp, "diagonalization")
     END IF
     IF (obj%conv_threshold_ispresent) THEN
        CALL xml_NewElement(xp, "conv_threshold")
           CALL xml_addCharacters(xp, obj%conv_threshold, fmt='s16')
        CALL xml_EndElement(xp, "conv_threshold")
     END IF
     IF (obj%q_gipaw_ispresent) THEN
        CALL xml_NewElement(xp, "q_gipaw")
           CALL xml_addCharacters(xp, obj%q_gipaw, fmt='s16')
        CALL xml_EndElement(xp, "q_gipaw")
     END IF
     IF (obj%r_rand_ispresent) THEN
        CALL xml_NewElement(xp, "r_rand")
           CALL xml_addCharacters(xp, obj%r_rand, fmt='s16')
        CALL xml_EndElement(xp, "r_rand")
     END IF
     IF (obj%spline_ps_ispresent) THEN
        CALL xml_NewElement(xp, "spline_ps")
           CALL xml_addCharacters(xp, obj%spline_ps)
        CALL xml_EndElement(xp, "spline_ps")
     END IF
     IF (obj%paw_proj_ispresent) THEN
        CALL xml_NewElement(xp, "paw_proj")
           CALL xml_addCharacters(xp, obj%paw_proj)
        CALL xml_EndElement(xp, "paw_proj")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_scf_gipaw

   SUBROUTINE qes_write_nmr(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(nmr_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%use_nmr_macroscopic_shape_ispresent) THEN
        CALL xml_NewElement(xp, "use_nmr_macroscopic_shape")
           CALL xml_addCharacters(xp, obj%use_nmr_macroscopic_shape)
        CALL xml_EndElement(xp, "use_nmr_macroscopic_shape")
     END IF
     IF (obj%nmr_macroscopic_shape_ispresent) THEN
        CALL xml_NewElement(xp, "nmr_macroscopic_shape")
           CALL xml_addCharacters(xp, obj%nmr_macroscopic_shape, fmt='s16')
        CALL xml_EndElement(xp, "nmr_macroscopic_shape")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_nmr

   SUBROUTINE qes_write_efg(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(efg_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%efg_q_ispresent) THEN
        CALL xml_NewElement(xp, "efg_q")
           CALL xml_addCharacters(xp, obj%efg_q, fmt='s16')
        CALL xml_EndElement(xp, "efg_q")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_efg

   SUBROUTINE qes_write_hfi(xp, obj) 
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(hfi_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%hfi_output_unit_ispresent) THEN
        CALL xml_NewElement(xp, "hfi_output_unit")
           CALL xml_addCharacters(xp, TRIM(obj%hfi_output_unit))
        CALL xml_EndElement(xp, "hfi_output_unit")
     END IF
     IF (obj%hfi_nuclear_g_tensor_ispresent) THEN
        CALL xml_NewElement(xp, "hfi_nuclear_g_tensor")
           CALL xml_addCharacters(xp, obj%hfi_nuclear_g_tensor, fmt='s16')
        CALL xml_EndElement(xp, "hfi_nuclear_g_tensor")
     END IF
     IF (obj%core_relax_method_ispresent) THEN
        CALL xml_NewElement(xp, "core_relax_method")
           CALL xml_addCharacters(xp, obj%core_relax_method)
        CALL xml_EndElement(xp, "core_relax_method")
     END IF
     IF (obj%hfi_via_reconstruction_only_ispresent) THEN
        CALL xml_NewElement(xp, "hfi_via_reconstruction_only")
           CALL xml_addCharacters(xp, obj%hfi_via_reconstruction_only)
        CALL xml_EndElement(xp, "hfi_via_reconstruction_only")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_hfi

  !
END MODULE qesgipaw_write_module