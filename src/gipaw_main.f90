!
! Copyright (C) 2001-2009 GIPAW and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM gipaw_main
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the magnetic response program. 
  ! ... It controls the initialization routines.
  ! ... Features: NMR chemical shifts
  !...            EPR g-tensor
  ! ... Ported to Espresso by:
  ! ... D. Ceresoli, A. P. Seitsonen, U. Gerstamnn and  F. Mauri
  ! ...
  ! ... References (NMR):
  ! ... F. Mauri and S. G. Louie Phys. Rev. Lett. 76, 4246 (1996)
  ! ... F. Mauri, B. G. Pfrommer, S. G. Louie, Phys. Rev. Lett. 77, 5300 (1996)
  ! ... T. Gregor, F. Mauri, and R. Car, J. Chem. Phys. 111, 1815 (1999)
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. B 63, 245101 (2001)
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. Lett. 91, 196401 (2003)
  ! ...
  ! ... References (g-tensor):
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. Lett. 88, 086403 (2002)
  ! ...
  USE kinds,           ONLY : DP
  USE io_files,        ONLY : tmp_dir, create_directory, prefix
  USE io_global,       ONLY : stdout, meta_ionode, meta_ionode_id
  USE mp,              ONLY : mp_bcast
  USE cell_base,       ONLY : tpiba
  USE gipaw_module,    ONLY : job, q_gipaw, max_seconds
  USE check_stop  ,    ONLY : check_stop_init
  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag
  USE mp_global,       ONLY : mp_startup, nproc_pool_file
  USE mp_world,        ONLY : world_comm
  USE mp_images,       ONLY : nimage, my_image_id
  USE mp_pools,        ONLY : intra_pool_comm
  USE mp_bands,        ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_bands,        ONLY : nbgrp
  USE mp_pools,        ONLY : nproc_pool
  USE environment,     ONLY : environment_start, environment_end
  USE lsda_mod,        ONLY : nspin
  USE wvfct,           ONLY : nbnd
  USE uspp,            ONLY : okvan
  USE wvfct,           ONLY : nbnd
  USE io_global,       ONLY : stdout
  USE noncollin_module,ONLY : noncolin
  USE cellmd,          ONLY : cell_factor
  USE xml_routines
  USE command_line_options, ONLY: input_file_, command_line, ndiag_
  ! for pluginization
  USE input_parameters, ONLY : nat_ => nat, ntyp_ => ntyp
  USE input_parameters, ONLY : assume_isolated_ => assume_isolated, &
                               ibrav_ => ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp
  USE cell_base,        ONLY : ibrav
  ! end
  USE gipaw_version
  !------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER (LEN=9)   :: code = 'GIPAW'
  CHARACTER (LEN=10)  :: dirname = 'dummy'
  LOGICAL, EXTERNAL  :: check_para_diag
  !------------------------------------------------------------------------

  ! begin with the initialization part
#ifdef __MPI
  call mp_startup(start_images=.true.)
#else
  call mp_startup(start_images=.false.)
#endif

!!#ifndef __BANDS
!!  call mp_start_diag(ndiag_, world_comm, intra_pool_comm, do_distr_diag_inside_bgrp_=.false.)
!!#else
!!  call mp_start_diag(ndiag_, world_comm, intra_pool_comm, do_distr_diag_inside_bgrp_=.true.)
!!#endif
  call set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm)

  call environment_start (code)

  ! read plugin command line arguments, if any
  if (meta_ionode) call plugin_arguments()
  call plugin_arguments_bcast( meta_ionode_id, world_comm )

#ifndef __BANDS
  if (nbgrp > 1) &
    call errore('gipaw_main', 'configure and recompile GIPAW with --enable-band-parallel', 1)
#endif

  write(stdout,*)
  write(stdout,'(5X,''***** This is GIPAW git revision '',A,'' *****'')') gipaw_git_revision
  write(stdout,'(5X,''***** you can cite: N. Varini et al., Comp. Phys. Comm. 184, 1827 (2013)  *****'')')
  write(stdout,'(5X,''***** in publications or presentations arising from this work.            *****'')')
  write(stdout,*)
 
  write(stdout,'(5X,''Parallelizing q-star over'',I2,'' images'')') nimage
  if (nimage > 7) write(stdout,'(5X,''ATTENTION: optimal number of images is 7'')')

  call gipaw_readin()
  call check_stop_init( max_seconds )

  io_level = 1
 
  ! read ground state wavefunctions
  cell_factor = 1.2
  call read_file

  call set_para_diag(nbnd, use_para_diag)
  
  call gipaw_openfil

  if (gamma_only) call errore ('gipaw_main', 'Cannot run GIPAW with gamma_only == .true. ', 1)
#ifdef __BANDS
  if (nbgrp > 1) &
    call errore('gipaw_main', 'Cannot use band-parallelization without wf_collect in SCF', 1)
#endif
  if (noncolin) call errore('gipaw_main', 'non-collinear not supported yet', 1)

  nat_ = nat
  ntyp_ = ntyp
  ibrav_ = ibrav
  assume_isolated_ = 'none'

  call plugin_read_input()
  call gipaw_allocate()
  call gipaw_setup()
  call gipaw_summary()
  
  ! convert q_gipaw into units of tpiba
  q_gipaw = q_gipaw / tpiba
  
  ! image initialization
  if (nimage > 1) then
     write(dirname, fmt='(I5.5)') my_image_id
     tmp_dir = trim(tmp_dir) // '/' // trim(dirname)
     call create_directory(tmp_dir)
  endif 

#ifdef __BANDS
  call init_parallel_over_band(inter_bgrp_comm, nbnd)
#endif

  ! calculation
  select case ( trim(job) )
  case ( 'nmr' )
     call suscept_crystal   
     
  case ( 'g_tensor' )
     if (nspin /= 2) call errore('gipaw_main', 'g_tensor is only for spin-polarized', 1)
     if (okvan) call errore('gipaw_main', 'g_tensor not available with ultrasoft', 1) 
     call suscept_crystal

  case ( 'knight' )
     call knight_shift
     
  case ( 'f-sum' )
     call suscept_crystal
     
  case ( 'efg' )
     call efg

  case ( 'mossbauer' )
     call mossbauer
     
  case ( 'hyperfine' )
     if (nspin /= 2) call errore('gipaw_main', 'hyperfine is only for spin-polarized', 1)
     call hyperfine
     
  case default
     call errore('gipaw_main', 'wrong or undefined job in input', 1)
  end select

  ! open XML output file
  call open_xml_output(trim(tmp_dir) // '/' // trim(prefix) // '-gipaw.xml')
  call xml_output_generalinfo
  call xml_output_parallelinfo
  call xml_output_namelist
  call xml_output_results
  call xml_output_status
  call xml_output_timinginfo
  call xml_output_closed
  call close_xml_output
  
  ! print timings and stop the code
  call gipaw_closefil
  call print_clock_gipaw
  call environment_end(code)
  call gipaw_deallocate
  call stop_code( .true. )
  
  STOP
  
END PROGRAM gipaw_main

