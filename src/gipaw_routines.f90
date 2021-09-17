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
  USE io_files,         ONLY : prefix, tmp_dir  
  USE io_global,        ONLY : ionode, stdout, qestdin
  USE uspp_data,        ONLY : spline_ps
  USE mp_images,        ONLY : my_image_id
  USE open_close_input_file

  ! -- local variables ---------------------------------------------------
  implicit none
  integer :: ios
  logical :: is_xml
  character(len=256), external :: trimcheck
  namelist /inputgipaw/ job, prefix, tmp_dir, conv_threshold, restart_mode, &
                        q_gipaw, iverbosity, filcurr, filfield, filnics, pawproj, &
                        use_nmr_macroscopic_shape, nmr_macroscopic_shape, &
                        spline_ps, isolve, q_efg, max_seconds, r_rand, &
                        hfi_output_unit, hfi_nuclear_g_factor, &
                        core_relax_method, diagonalization, verbosity, &
                        hfi_via_reconstruction_only

  if (.not. ionode .or. my_image_id > 0) goto 400
    
  !!call input_from_file()
  ios = open_input_file(is_xml=is_xml)
  if (ios > 0) call errore('gipaw_readin', 'error reading input file', 1)
  if (is_xml) call errore('gipaw_readin', 'XML not yet supported', 1)

  ! define input default values
  call get_environment_variable( 'ESPRESSO_TMPDIR', tmp_dir ) 
  if (trim(tmp_dir) == ' ') tmp_dir = './scratch/'
  tmp_dir = trimcheck(tmp_dir)
  job = ''
  prefix = 'pwscf'
  restart_mode = 'restart'
  conv_threshold = 1d-14
  q_gipaw = 0.01d0
  iverbosity = -1
  verbosity = 'low'
  filcurr = ''
  filfield = ''
  filnics = ''
  use_nmr_macroscopic_shape = .true.
  nmr_macroscopic_shape(:,:) = 0.d0
  nmr_macroscopic_shape(1,1) = 2.d0 / 3.d0
  nmr_macroscopic_shape(2,2) = 2.d0 / 3.d0
  nmr_macroscopic_shape(3,3) = 2.d0 / 3.d0
  spline_ps = .true.
  isolve = -1
  diagonalization = 'david'
  core_relax_method = 1
  hfi_via_reconstruction_only = .false.

  hfi_output_unit = 'MHz'
  hfi_nuclear_g_factor(:) = 1.0
  q_efg(:)  = 1.0
  pawproj(:) = .false.
  r_rand = 0.1
  max_seconds  =  1.d7

  print*, qestdin
  read(qestdin, inputgipaw, err=200, iostat=ios)

  ! check input
  if (max_seconds < 0.1d0) call errore ('gipaw_readin', ' wrong max_seconds', 1)
200 call errore('gipaw_readin', 'reading inputgipaw namelist', abs(ios))

  ! further checks
  if (isolve /= -1) &
     call errore('gipaw_readin', '*** isolve is obsolete, use diagonalization instead ***', 1)
  if (iverbosity /= -1) &
     call errore('gipaw_readin', '*** iverbosity is obsolete, use verbosity instead ***', 1)

  select case (diagonalization)
     case('david')
       isolve = 0
     case('cg')
       isolve = 1
     !!case('ppcg')
     !!  isolve = 2
     !!case('rmm-diis')
     !!  isolve = 3
     case default
       call errore('gipaw_readin', 'diagonalization can be ''david'', ''cg'' or ''ppcg''', 1)
  end select

  select case (verbosity)
     case('low')
       iverbosity = 1
     case('medium')
       iverbosity = 11
     case('high')
       iverbosity = 21
     case default
       call errore('gipaw_readin', 'verbosity can be ''low'', ''medium'' or ''high''', 1)
  end select

  ! dump input file in output
  if (verbosity == 'high') call dump_input_file

  ! close input file
  ios = close_input_file()

400 continue

#ifdef __MPI
  ! broadcast input variables  
  call gipaw_bcast_input
#endif

END SUBROUTINE gipaw_readin


!-----------------------------------------------------------------------
SUBROUTINE dump_input_file
  !-----------------------------------------------------------------------
  !
  ! ... dump input file to the output
  !
  USE io_global,        ONLY : ionode, stdout, qestdin
  USE gipaw_module
  implicit none
  character(len=256) :: line
  integer :: ios

  if (.not. ionode) return

  write(stdout,*)
  write(stdout,'(5X,''------------------- Input file: --------------------'')')
  rewind(qestdin)
  do
      read(qestdin,'(A)',iostat=ios) line
      if (ios == 0) then
          write(stdout,'(A)') trim(line)
      else
          exit
      endif
  end do
  write(stdout,*)

END SUBROUTINE dump_input_file


#ifdef __MPI
!-----------------------------------------------------------------------
SUBROUTINE gipaw_bcast_input
  !-----------------------------------------------------------------------
  !
  ! ... Broadcast input data to all processors 
  !
  USE gipaw_module
  USE mp_world,      ONLY : world_comm
  USE mp,            ONLY : mp_bcast
  USE io_files,      ONLY : prefix, tmp_dir
  USE uspp_data,     ONLY : spline_ps

  implicit none
  integer :: root = 0

  call mp_bcast(job, root, world_comm)
  call mp_bcast(prefix, root, world_comm)
  call mp_bcast(tmp_dir, root, world_comm)
  call mp_bcast(conv_threshold, root, world_comm)
  call mp_bcast(q_gipaw, root, world_comm)
  call mp_bcast(iverbosity, root, world_comm)
  call mp_bcast(filcurr, root, world_comm)
  call mp_bcast(filfield, root, world_comm)
  call mp_bcast(filnics, root, world_comm)
  call mp_bcast(use_nmr_macroscopic_shape, root, world_comm)
  call mp_bcast(nmr_macroscopic_shape, root, world_comm)
  call mp_bcast(spline_ps, root, world_comm)
  call mp_bcast(isolve, root, world_comm)
  call mp_bcast(core_relax_method, root, world_comm)
  call mp_bcast(hfi_via_reconstruction_only, root, world_comm)
  call mp_bcast(hfi_output_unit, root, world_comm)
  call mp_bcast(hfi_nuclear_g_factor, root, world_comm)
  call mp_bcast(q_efg, root, world_comm)
  call mp_bcast(pawproj, root, world_comm)
  call mp_bcast(r_rand, root, world_comm)
  call mp_bcast(max_seconds, root, world_comm)
  call mp_bcast(restart_mode, root, world_comm)
  call mp_bcast(verbosity, root, world_comm)
  call mp_bcast(diagonalization, root, world_comm)

END SUBROUTINE gipaw_bcast_input
#endif



!-----------------------------------------------------------------------
SUBROUTINE gipaw_allocate
  !-----------------------------------------------------------------------
  !
  ! ... Allocate memory for GIPAW
  !
  USE gipaw_module
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
  USE cellmd,        ONLY : cell_factor
  USE gvecw,         ONLY : ecutwfc
  USE uspp_data,     ONLY : spline_ps
  implicit none

  if (.not. spline_ps) then
      write(stdout,*)
      call infomsg('gipaw_summary', 'spline_ps is .false., expect some extrapolation errors')
  endif

  write(stdout,*)
  write(stdout,"(5X,'q-space interpolation up to ',F8.2,' Rydberg')") ecutwfc*cell_factor
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

  call plugin_summary

  flush(stdout)

END SUBROUTINE gipaw_summary
  

!-----------------------------------------------------------------------
SUBROUTINE gipaw_openfil
  !-----------------------------------------------------------------------
  !
  ! ... Open files needed for GIPAW
  !
  USE gipaw_module
  USE wvfct,            ONLY : nbnd, npwx
  USE io_files,         ONLY : iunwfc, nwordwfc
  USE noncollin_module, ONLY : npol
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

END SUBROUTINE gipaw_openfil


!-----------------------------------------------------------------------
SUBROUTINE gipaw_closefil
  !-----------------------------------------------------------------------
  !
  ! ... Close files opened by GIPAW, if any
  !
  return

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

  write(stdout,*) '    Plugins'
  call plugin_clock()
  write(stdout,*)

END SUBROUTINE print_clock_gipaw


!-----------------------------------------------------------------------
SUBROUTINE gipaw_memory_report
  !-----------------------------------------------------------------------
  !
  ! ... Print estimated memory usage
  !
  USE io_global,                 ONLY : stdout
  USE noncollin_module,          ONLY : npol
  USE uspp,                      ONLY : okvan, nkb
  USE fft_base,                  ONLY : dffts
  USE paw_gipaw,                 ONLY : paw_nkb
  USE pwcom
  IMPLICIT NONE
  integer, parameter :: Mb=1024*1024, complex_size=16, real_size=8

  ! the conversions to double prevent integer overflow in very large run
  write(stdout,'(5x,"Largest allocated arrays",5x,"est. size (Mb)",5x,"dimensions")')

  write(stdout,'(8x,"KS wavefunctions at k     ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     complex_size*nbnd*npol*DBLE(npwx)/Mb, npwx*npol,nbnd
  write(stdout,'(8x,"KS wavefunctions at k+q   ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     complex_size*nbnd*npol*DBLE(npwx)/Mb, npwx*npol,nbnd
  write(stdout,'(8x,"First-order wavefunctions ",f10.2," Mb",5x,"(",i8,",",i5,",",i3")")') &
     complex_size*nbnd*npol*DBLE(npwx)*10/Mb, npwx*npol,nbnd,10
  if (okvan) &
  write(stdout,'(8x,"First-order wavefunct (US)",f10.2," Mb",5x,"(",i8,",",i5,",",i3")")') &
     complex_size*nbnd*npol*DBLE(npwx)*6/Mb, npwx*npol,nbnd,6

  write(stdout,'(8x,"Charge/spin density       ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     real_size*dble(dffts%nnr)*nspin/Mb, dffts%nnr, nspin
  write(stdout,'(8x,"Induced current           ",f10.2," Mb",5x,"(",i8,",",i5,",",i1,",",i1")")') &
     real_size*dble(dffts%nnr)*9*nspin/Mb, dffts%nnr,3,3,nspin
  write(stdout,'(8x,"Induced magnetic field    ",f10.2," Mb",5x,"(",i8,",",i5,",",i1,",",i1")")') &
     real_size*dble(dffts%nnr)*9*nspin/Mb, dffts%nnr,3,3,nspin
  
  write(stdout,'(8x,"NL pseudopotentials       ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     complex_size*nkb*DBLE(npwx)/Mb, npwx, nkb
  write(stdout,'(8x,"GIPAW NL terms            ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     complex_size*paw_nkb*DBLE(npwx)/Mb, npwx, paw_nkb
  write(stdout,*)

END SUBROUTINE gipaw_memory_report


