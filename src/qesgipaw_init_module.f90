!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE qesgipaw_init_module
  !
  ! Auto-generated code: don't edit or at least don't commit changes
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes_1.0
  !
  USE kinds, only: DP
  USE qesgipaw_types_module
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: qesgipaw_init
  !
  INTERFACE qesgipaw_init
    !
    MODULE PROCEDURE qes_init_gipaw
    MODULE PROCEDURE qes_init_general_info
    MODULE PROCEDURE qes_init_format
    MODULE PROCEDURE qes_init_creator
    MODULE PROCEDURE qes_init_created
    MODULE PROCEDURE qes_init_parallel_info
    MODULE PROCEDURE qes_init_closed
    MODULE PROCEDURE qes_init_inputGIPAW
    MODULE PROCEDURE qes_init_job
    MODULE PROCEDURE qes_init_files
    MODULE PROCEDURE qes_init_scf_gipaw
    MODULE PROCEDURE qes_init_nmr
    MODULE PROCEDURE qes_init_efg
    MODULE PROCEDURE qes_init_hfi
    !
  END INTERFACE qesgipaw_init
  !
  CONTAINS
  !
  !
  SUBROUTINE qes_init_gipaw(obj, tagname, inputGIPAW, general_info, parallel_info, status, cputime, closed)
    !
    IMPLICIT NONE
    !
    TYPE(gipaw_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(general_info_type),OPTIONAL,INTENT(IN) :: general_info
    TYPE(parallel_info_type),OPTIONAL,INTENT(IN) :: parallel_info
    TYPE(inputGIPAW_type),INTENT(IN) :: inputGIPAW
    INTEGER,OPTIONAL,INTENT(IN) :: status
    INTEGER,OPTIONAL,INTENT(IN) :: cputime
    TYPE(closed_type),OPTIONAL,INTENT(IN) :: closed
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(general_info)) THEN 
      obj%general_info_ispresent = .TRUE. 
      obj%general_info = general_info
    ELSE 
      obj%general_info_ispresent = .FALSE.
    END IF
    IF ( PRESENT(parallel_info)) THEN 
      obj%parallel_info_ispresent = .TRUE. 
      obj%parallel_info = parallel_info
    ELSE 
      obj%parallel_info_ispresent = .FALSE.
    END IF
    obj%inputGIPAW = inputGIPAW
    IF ( PRESENT(status)) THEN 
      obj%status_ispresent = .TRUE. 
      obj%status = status
    ELSE 
      obj%status_ispresent = .FALSE.
    END IF
    IF ( PRESENT(cputime)) THEN 
      obj%cputime_ispresent = .TRUE. 
      obj%cputime = cputime
    ELSE 
      obj%cputime_ispresent = .FALSE.
    END IF
    IF ( PRESENT(closed)) THEN 
      obj%closed_ispresent = .TRUE. 
      obj%closed = closed
    ELSE 
      obj%closed_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_gipaw 
  !
  !
  SUBROUTINE qes_init_general_info(obj, tagname, format, creator, created, job)
    !
    IMPLICIT NONE
    !
    TYPE(general_info_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(format_type),INTENT(IN) :: format
    TYPE(creator_type),INTENT(IN) :: creator
    TYPE(created_type),INTENT(IN) :: created
    CHARACTER(LEN=*),INTENT(IN) :: job
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%format = format
    obj%creator = creator
    obj%created = created
    obj%job = job
    !
  END SUBROUTINE qes_init_general_info 
  !
  !
  SUBROUTINE qes_init_format(obj, tagname, NAME, VERSION, format)
    !
    IMPLICIT NONE
    !
    TYPE(format_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    CHARACTER(LEN=*), INTENT(IN) :: VERSION
    CHARACTER(LEN=*), INTENT(IN) :: format
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%NAME = NAME
    obj%VERSION = VERSION
    !
    obj%format = format
    !
  END SUBROUTINE qes_init_format 
  !
  !
  SUBROUTINE qes_init_creator(obj, tagname, NAME, VERSION, creator)
    !
    IMPLICIT NONE
    !
    TYPE(creator_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    CHARACTER(LEN=*), INTENT(IN) :: VERSION
    CHARACTER(LEN=*), INTENT(IN) :: creator
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%NAME = NAME
    obj%VERSION = VERSION
    !
    obj%creator = creator
    !
  END SUBROUTINE qes_init_creator 
  !
  !
  SUBROUTINE qes_init_created(obj, tagname, DATE, TIME, created)
    !
    IMPLICIT NONE
    !
    TYPE(created_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: DATE
    CHARACTER(LEN=*), INTENT(IN) :: TIME
    CHARACTER(LEN=*), INTENT(IN) :: created
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%DATE = DATE
    obj%TIME = TIME
    !
    obj%created = created
    !
  END SUBROUTINE qes_init_created 
  !
  !
  SUBROUTINE qes_init_parallel_info(obj, tagname, nprocs, nthreads, ntasks, nbgrp, npool, ndiag, nimages)
    !
    IMPLICIT NONE
    !
    TYPE(parallel_info_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,INTENT(IN) :: nprocs
    INTEGER,INTENT(IN) :: nthreads
    INTEGER,INTENT(IN) :: ntasks
    INTEGER,INTENT(IN) :: nbgrp
    INTEGER,INTENT(IN) :: npool
    INTEGER,INTENT(IN) :: ndiag
    INTEGER,OPTIONAL,INTENT(IN) :: nimages
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%nprocs = nprocs
    obj%nthreads = nthreads
    obj%ntasks = ntasks
    obj%nbgrp = nbgrp
    obj%npool = npool
    obj%ndiag = ndiag
    IF ( PRESENT(nimages)) THEN 
      obj%nimages_ispresent = .TRUE. 
      obj%nimages = nimages
    ELSE 
      obj%nimages_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_parallel_info 
  !
  !
  SUBROUTINE qes_init_closed(obj, tagname, DATE, TIME, closed)
    !
    IMPLICIT NONE
    !
    TYPE(closed_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: DATE
    CHARACTER(LEN=*), INTENT(IN) :: TIME
    CHARACTER(LEN=*), INTENT(IN) :: closed
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%DATE = DATE
    obj%TIME = TIME
    !
    obj%closed = closed
    !
  END SUBROUTINE qes_init_closed 
  !
  !
  SUBROUTINE qes_init_inputGIPAW(obj, tagname, control_job, scf_gipaw, files, nmr_input, efg_input, hfi_input)
    !
    IMPLICIT NONE
    !
    TYPE(inputGIPAW_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(job_type),OPTIONAL,INTENT(IN) :: control_job
    TYPE(scf_gipaw_type),OPTIONAL,INTENT(IN) :: scf_gipaw
    TYPE(files_type),OPTIONAL,INTENT(IN) :: files
    TYPE(nmr_type),OPTIONAL,INTENT(IN) :: nmr_input
    TYPE(efg_type),OPTIONAL,INTENT(IN) :: efg_input
    TYPE(hfi_type),OPTIONAL,INTENT(IN) :: hfi_input
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(control_job)) THEN 
      obj%control_job_ispresent = .TRUE. 
      obj%control_job = control_job
    ELSE 
      obj%control_job_ispresent = .FALSE.
    END IF
    IF ( PRESENT(scf_gipaw)) THEN 
      obj%scf_gipaw_ispresent = .TRUE. 
      obj%scf_gipaw = scf_gipaw
    ELSE 
      obj%scf_gipaw_ispresent = .FALSE.
    END IF
    IF ( PRESENT(files)) THEN 
      obj%files_ispresent = .TRUE. 
      obj%files = files
    ELSE 
      obj%files_ispresent = .FALSE.
    END IF
    IF ( PRESENT(nmr_input)) THEN 
      obj%nmr_input_ispresent = .TRUE. 
      obj%nmr_input = nmr_input
    ELSE 
      obj%nmr_input_ispresent = .FALSE.
    END IF
    IF ( PRESENT(efg_input)) THEN 
      obj%efg_input_ispresent = .TRUE. 
      obj%efg_input = efg_input
    ELSE 
      obj%efg_input_ispresent = .FALSE.
    END IF
    IF ( PRESENT(hfi_input)) THEN 
      obj%hfi_input_ispresent = .TRUE. 
      obj%hfi_input = hfi_input
    ELSE 
      obj%hfi_input_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_inputGIPAW 
  !
  !
  SUBROUTINE qes_init_job(obj, tagname, job, restart_mode, verbosity, max_seconds)
    !
    IMPLICIT NONE
    !
    TYPE(job_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: job
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: restart_mode
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: verbosity
    REAL(DP),OPTIONAL,INTENT(IN) :: max_seconds
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(job)) THEN 
      obj%job_ispresent = .TRUE. 
      obj%job = job
    ELSE 
      obj%job_ispresent = .FALSE.
    END IF
    IF ( PRESENT(restart_mode)) THEN 
      obj%restart_mode_ispresent = .TRUE. 
      obj%restart_mode = restart_mode
    ELSE 
      obj%restart_mode_ispresent = .FALSE.
    END IF
    IF ( PRESENT(verbosity)) THEN 
      obj%verbosity_ispresent = .TRUE. 
      obj%verbosity = verbosity
    ELSE 
      obj%verbosity_ispresent = .FALSE.
    END IF
    IF ( PRESENT(max_seconds)) THEN 
      obj%max_seconds_ispresent = .TRUE. 
      obj%max_seconds = max_seconds
    ELSE 
      obj%max_seconds_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_job 
  !
  !
  SUBROUTINE qes_init_files(obj, tagname, prefix, tmp_dir, filcurr, fildfield, filnics)
    !
    IMPLICIT NONE
    !
    TYPE(files_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: prefix
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: tmp_dir
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filcurr
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: fildfield
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: filnics
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(prefix)) THEN 
      obj%prefix_ispresent = .TRUE. 
      obj%prefix = prefix
    ELSE 
      obj%prefix_ispresent = .FALSE.
    END IF
    IF ( PRESENT(tmp_dir)) THEN 
      obj%tmp_dir_ispresent = .TRUE. 
      obj%tmp_dir = tmp_dir
    ELSE 
      obj%tmp_dir_ispresent = .FALSE.
    END IF
    IF ( PRESENT(filcurr)) THEN 
      obj%filcurr_ispresent = .TRUE. 
      obj%filcurr = filcurr
    ELSE 
      obj%filcurr_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fildfield)) THEN 
      obj%fildfield_ispresent = .TRUE. 
      obj%fildfield = fildfield
    ELSE 
      obj%fildfield_ispresent = .FALSE.
    END IF
    IF ( PRESENT(filnics)) THEN 
      obj%filnics_ispresent = .TRUE. 
      obj%filnics = filnics
    ELSE 
      obj%filnics_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_files 
  !
  !
  SUBROUTINE qes_init_scf_gipaw(obj, tagname, diagonalization, conv_threshold, q_gipaw, r_rand,&
                               spline_ps, paw_proj)
    !
    IMPLICIT NONE
    !
    TYPE(scf_gipaw_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: diagonalization
    REAL(DP),OPTIONAL,INTENT(IN) :: conv_threshold
    REAL(DP),OPTIONAL,INTENT(IN) :: q_gipaw
    REAL(DP),OPTIONAL,INTENT(IN) :: r_rand
    LOGICAL,OPTIONAL,INTENT(IN) :: spline_ps
    LOGICAL, DIMENSION(:),OPTIONAL,INTENT(IN) :: paw_proj
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(diagonalization)) THEN 
      obj%diagonalization_ispresent = .TRUE. 
      obj%diagonalization = diagonalization
    ELSE 
      obj%diagonalization_ispresent = .FALSE.
    END IF
    IF ( PRESENT(conv_threshold)) THEN 
      obj%conv_threshold_ispresent = .TRUE. 
      obj%conv_threshold = conv_threshold
    ELSE 
      obj%conv_threshold_ispresent = .FALSE.
    END IF
    IF ( PRESENT(q_gipaw)) THEN 
      obj%q_gipaw_ispresent = .TRUE. 
      obj%q_gipaw = q_gipaw
    ELSE 
      obj%q_gipaw_ispresent = .FALSE.
    END IF
    IF ( PRESENT(r_rand)) THEN 
      obj%r_rand_ispresent = .TRUE. 
      obj%r_rand = r_rand
    ELSE 
      obj%r_rand_ispresent = .FALSE.
    END IF
    IF ( PRESENT(spline_ps)) THEN 
      obj%spline_ps_ispresent = .TRUE. 
      obj%spline_ps = spline_ps
    ELSE 
      obj%spline_ps_ispresent = .FALSE.
    END IF
    IF ( PRESENT(paw_proj)) THEN 
      obj%paw_proj_ispresent = .TRUE. 
      obj%paw_proj = paw_proj
    ELSE 
      obj%paw_proj_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_scf_gipaw 
  !
  !
  SUBROUTINE qes_init_nmr(obj, tagname, use_nmr_macroscopic_shape, nmr_macroscopic_shape)
    !
    IMPLICIT NONE
    !
    TYPE(nmr_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,OPTIONAL,INTENT(IN) :: use_nmr_macroscopic_shape
    REAL(DP), DIMENSION(9),OPTIONAL,INTENT(IN) :: nmr_macroscopic_shape
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(use_nmr_macroscopic_shape)) THEN 
      obj%use_nmr_macroscopic_shape_ispresent = .TRUE. 
      obj%use_nmr_macroscopic_shape = use_nmr_macroscopic_shape
    ELSE 
      obj%use_nmr_macroscopic_shape_ispresent = .FALSE.
    END IF
    IF ( PRESENT(nmr_macroscopic_shape)) THEN 
      obj%nmr_macroscopic_shape_ispresent = .TRUE. 
      obj%nmr_macroscopic_shape = nmr_macroscopic_shape
    ELSE 
      obj%nmr_macroscopic_shape_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_nmr 
  !
  !
  SUBROUTINE qes_init_efg(obj, tagname, efg_q)
    !
    IMPLICIT NONE
    !
    TYPE(efg_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), DIMENSION(:),OPTIONAL,INTENT(IN) :: efg_q
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(efg_q)) THEN 
      obj%efg_q_ispresent = .TRUE. 
      obj%efg_q = efg_q
    ELSE 
      obj%efg_q_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_efg 
  !
  !
  SUBROUTINE qes_init_hfi(obj, tagname, hfi_output_unit, hfi_nuclear_g_tensor, core_relax_method, hfi_via_reconstruction_only)
    !
    IMPLICIT NONE
    !
    TYPE(hfi_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: hfi_output_unit
    REAL(DP), DIMENSION(:),OPTIONAL,INTENT(IN) :: hfi_nuclear_g_tensor
    INTEGER,OPTIONAL,INTENT(IN) :: core_relax_method
    LOGICAL,OPTIONAL,INTENT(IN) :: hfi_via_reconstruction_only
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(hfi_output_unit)) THEN 
      obj%hfi_output_unit_ispresent = .TRUE. 
      obj%hfi_output_unit = hfi_output_unit
    ELSE 
      obj%hfi_output_unit_ispresent = .FALSE.
    END IF
    IF ( PRESENT(hfi_nuclear_g_tensor)) THEN 
      obj%hfi_nuclear_g_tensor_ispresent = .TRUE. 
      obj%hfi_nuclear_g_tensor = hfi_nuclear_g_tensor
    ELSE 
      obj%hfi_nuclear_g_tensor_ispresent = .FALSE.
    END IF
    IF ( PRESENT(core_relax_method)) THEN 
      obj%core_relax_method_ispresent = .TRUE. 
      obj%core_relax_method = core_relax_method
    ELSE 
      obj%core_relax_method_ispresent = .FALSE.
    END IF
    IF ( PRESENT(hfi_via_reconstruction_only)) THEN 
      obj%hfi_via_reconstruction_only_ispresent = .TRUE. 
      obj%hfi_via_reconstruction_only = hfi_via_reconstruction_only
    ELSE 
      obj%hfi_via_reconstruction_only_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_hfi 
  !
  !
END MODULE qesgipaw_init_module