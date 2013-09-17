!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE gipaw_readin()
  !-----------------------------------------------------------------------
  !
  ! ... Read in the gipaw input file. The input file consists of a
  ! ... single namelist &inputgipaw. See doc/user-manual.pdf for the
  ! ... list of input keywords.
  !
  USE gipaw_module
  USE io_files,         ONLY : nd_nmbr, prefix, tmp_dir  
  USE io_global,        ONLY : ionode
  USE us,               ONLY : spline_ps
  USE input_parameters, ONLY : max_seconds
  USE mp_images,        ONLY : my_image_id

  ! -- local variables ---------------------------------------------------
  implicit none
  integer :: ios
  character(len=256), external :: trimcheck
  namelist /inputgipaw/ job, prefix, tmp_dir, conv_threshold, restart_mode, &
                        q_gipaw, iverbosity, filcurr, filfield, filnics, pawproj, &
                        use_nmr_macroscopic_shape, nmr_macroscopic_shape, &
                        spline_ps, isolve, q_efg, max_seconds, r_rand, &
                        hfi_output_unit, hfi_nuclear_g_factor, &
                        core_relax_method

  if (.not. ionode .or. my_image_id > 0) goto 400
    
  call input_from_file()
    
  ! define input default values
  call get_env( 'ESPRESSO_TMPDIR', tmp_dir ) 
  if (trim(tmp_dir) == ' ') tmp_dir = './scratch/'
  tmp_dir = trimcheck(tmp_dir)
  job = ''
  prefix = 'pwscf'
  restart_mode = 'restart'
  conv_threshold = 1d-14
  q_gipaw = 0.01d0
  iverbosity = 0
  filcurr = ''
  filfield = ''
  filnics = ''
  use_nmr_macroscopic_shape = .true.
  nmr_macroscopic_shape(:,:) = 0.d0
  nmr_macroscopic_shape(1,1) = 2.d0 / 3.d0
  nmr_macroscopic_shape(2,2) = 2.d0 / 3.d0
  nmr_macroscopic_shape(3,3) = 2.d0 / 3.d0
  spline_ps = .true.
  isolve = 0
  core_relax_method = 1

  hfi_output_unit = 'MHz'
  hfi_nuclear_g_factor(:) = 1.0
  q_efg(:)  = 1.0
  pawproj(:) = .false.
  r_rand = 0.1
  max_seconds  =  1.d7

  ! read input    
  read( 5, inputgipaw, err = 200, iostat = ios )
    
  ! check input
  if (max_seconds < 0.1d0) call errore ('gipaw_readin', ' wrong max_seconds', 1)
200 call errore( 'gipaw_readin', 'reading inputgipaw namelist', abs( ios ) )
    
400 continue
 
#ifdef __MPI
  ! broadcast input variables  
  call gipaw_bcast_input
#endif
    
END SUBROUTINE gipaw_readin


#ifdef __MPI
!-----------------------------------------------------------------------
SUBROUTINE gipaw_bcast_input
  !-----------------------------------------------------------------------
  !
  ! ... Broadcast input data to all processors 
  !
  USE gipaw_module
  USE mp,            ONLY : mp_bcast
  USE io_files,      ONLY : prefix, tmp_dir
  USE us,            ONLY : spline_ps
  USE input_parameters, ONLY : max_seconds

  implicit none
  integer :: root = 0

  call mp_bcast(job, root)
  call mp_bcast(prefix, root)
  call mp_bcast(tmp_dir, root)
  call mp_bcast(conv_threshold, root)
  call mp_bcast(q_gipaw, root)
  call mp_bcast(iverbosity, root)
  call mp_bcast(filcurr, root)
  call mp_bcast(filfield, root)
  call mp_bcast(filnics, root)
  call mp_bcast(use_nmr_macroscopic_shape, root)
  call mp_bcast(nmr_macroscopic_shape, root)
  call mp_bcast(spline_ps, root)
  call mp_bcast(isolve, root)
  call mp_bcast(core_relax_method, root)
  call mp_bcast(hfi_output_unit, root)
  call mp_bcast(hfi_nuclear_g_factor, root)
  call mp_bcast(q_efg, root)
  call mp_bcast(pawproj, root)
  call mp_bcast(r_rand, root)
  call mp_bcast(max_seconds, root)
  call mp_bcast(restart_mode, root)

END SUBROUTINE gipaw_bcast_input
#endif



!-----------------------------------------------------------------------
SUBROUTINE gipaw_allocate
  !-----------------------------------------------------------------------
  !
  ! ... Allocate memory for GIPAW
  !
  USE gipaw_module
  USE lsda_mod,      ONLY : nspin, lsda
  USE ions_base,     ONLY : ntyp => nsp
  USE paw_gipaw,     ONLY : paw_recon
  USE pwcom
    
  implicit none
  
  ! wavefunction at k+q  
  allocate(evq(npwx,nbnd))

  ! eigenvalues
  allocate(etq(nbnd,nkstot))

  ! GIPAW projectors
  if (.not. allocated(paw_recon)) allocate(paw_recon(ntyp))
    
END SUBROUTINE gipaw_allocate
  


!-----------------------------------------------------------------------
SUBROUTINE gipaw_summary
  !-----------------------------------------------------------------------
  !
  ! ... Print a short summary of the calculation
  !
  USE gipaw_module
  USE io_global,     ONLY : stdout
  USE ldaU,          ONLY : lda_plus_U
  implicit none

  write(stdout,*)

  write(stdout,'(5X,''GIPAW job: '',A)') job
  if (job(1:3) == 'nmr') then
    if (use_nmr_macroscopic_shape) then
      write(stdout,'(5X,''NMR macroscopic correction: yes'')')
      write(stdout,tens_fmt) nmr_macroscopic_shape
    else
      write(stdout,'(5X,''NMR macroscopic correction: no'')')
    endif
  endif

  if (lda_plus_U .eqv. .true.) call write_ns

  write(stdout,*)

  call flush_unit( stdout )

END SUBROUTINE gipaw_summary
  


!-----------------------------------------------------------------------
SUBROUTINE gipaw_openfil
  !-----------------------------------------------------------------------
  !
  ! ... Open files needed for GIPAW
  !
  USE gipaw_module
  USE io_global,        ONLY : stdout
  USE wvfct,            ONLY : nbnd, npwx
  USE ldaU,             ONLY : lda_plus_U, nwfcU
  USE klist,            ONLY : nks
  USE io_files,         ONLY : prefix, iunhub, iunwfc, &
                               nwordwfcU, nwordwfc, nwordatwfc, seqopn
  USE noncollin_module, ONLY : npol
  USE mp_global,        ONLY : kunit
  USE buffers,          ONLY : open_buffer
  USE control_flags,    ONLY : io_level    
  IMPLICIT NONE  

  logical :: exst

  !
  ! ... nwordwfc is the record length (IN REAL WORDS)
  ! ... for the direct-access file containing wavefunctions
  ! ... io_level > 0 : open a file; io_level <= 0 : open a buffer
  !
  nwordwfc = nbnd*npwx*npol
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst )

  ! ... Needed for LDA+U
  ! ... iunhub contains the (orthogonalized) atomic wfcs * S
  
  nwordwfcU = npwx*nwfcU*npol
  IF ( lda_plus_u ) &
     CALL open_buffer( iunhub, 'hub', nwordwfcU, io_level, exst )

END SUBROUTINE gipaw_openfil

!-----------------------------------------------------------------------
SUBROUTINE gipaw_closefil
  !-----------------------------------------------------------------------
  !
  ! ... Close files opened by GIPAW
  !
  USE ldaU,             ONLY : lda_plus_U  
  USE io_files,         ONLY : prefix, iunhub, iunwfc
  USE buffers,          ONLY : close_buffer

  call close_buffer( iunwfc, 'keep' )
  if ( lda_plus_u ) call close_buffer ( iunhub, status = 'keep' )

END SUBROUTINE gipaw_closefil



!-----------------------------------------------------------------------
SUBROUTINE print_clock_gipaw
  !-----------------------------------------------------------------------
  !
  ! ... Print clocks
  !
  USE io_global,  ONLY : stdout
  IMPLICIT NONE

  write(stdout,*) '    Initialization:'
  call print_clock ('gipaw_setup')
  write(stdout,*)
  write(stdout,*) '    Linear response'
  call print_clock ('greenf')
  call print_clock ('cgsolve')
  call print_clock ('ch_psi')
  call print_clock ('h_psiq')
  call print_clock ('u_kq')
  write(stdout,*)
  write(stdout,*) '    Apply operators'
  call print_clock ('h_psi')
  call print_clock ('apply_vel')
  write(stdout,*)
  write(stdout,*) '    Induced current'
  call print_clock ('j_para')
  call print_clock ('biot_savart')
  call print_clock ('c_sigma')
  write(stdout,*)
  write(stdout,*) '    Other routines'
  call print_clock ('f-sum')
  call print_clock ('efg')
  call print_clock ('hyperfine')
  call print_clock ('core_relax')
  write(stdout,*)
  write(stdout,*) '    General routines'
  call print_clock ('calbec')
  call print_clock ('fft')
  call print_clock ('ffts')
  call print_clock ('fftw')
  call print_clock ('cinterpolate')
  call print_clock ('davcio')
  call print_clock ('write_rec')
  write(stdout,*)

#ifdef __MPI
  write(stdout,*) '    Parallel routines'
  call print_clock ('reduce')  
  call print_clock( 'fft_scatter' )
  call print_clock( 'ALLTOALL' )
  write(stdout,*)
#endif
  call print_clock ('GIPAW') 

END SUBROUTINE print_clock_gipaw

