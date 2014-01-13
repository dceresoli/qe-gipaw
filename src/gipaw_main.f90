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
  USE io_files,        ONLY : prefix, tmp_dir
  USE io_global,       ONLY : stdout
  USE klist,           ONLY : nks
  USE mp,              ONLY : mp_bcast
  USE cell_base,       ONLY : tpiba
  USE cellmd,          ONLY : cell_factor
  USE gipaw_module,    ONLY : job, q_gipaw
  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag, twfcollect
  USE mp_global,       ONLY : mp_startup, nproc_pool_file
  USE mp_images,       ONLY : nimage, my_image_id
  USE mp_bands,        ONLY : inter_bgrp_comm, nbgrp
  USE mp_pools,        ONLY : nproc_pool
  USE check_stop,      ONLY : check_stop_init
  USE environment,     ONLY : environment_start
  USE lsda_mod,        ONLY : nspin
  USE wvfct,           ONLY : nbnd
  USE uspp,            ONLY : okvan
  USE wvfct,           ONLY : nbnd, npw 
  USE io_global,       ONLY : stdout
  USE fft_base,        ONLY : dffts
  USE noncollin_module,ONLY : noncolin
  USE gipaw_version
  USE iotk_module  
  USE xml_io_base
  !------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER (LEN=9)   :: code = 'QE'
  CHARACTER (LEN=10)  :: dirname = 'dummy'
  LOGICAL, EXTERNAL  :: check_para_diag
  !------------------------------------------------------------------------

  ! begin with the initialization part
#ifdef __MPI
  call mp_startup(start_images=.true.)
#else
  call mp_startup(start_images=.false.)
#endif
  call environment_start (code)

#ifndef __BANDS
  if (nbgrp > 1) &
    call errore('gipaw_main', 'configure and recompile GIPAW with --enable-band-parallel', 1)
#endif

  write(stdout,*)
  write(stdout,'(5X,''***** This is GIPAW svn revision '',A,'' *****'')') gipaw_svn_revision
  write(stdout,*)
 
  write(stdout,'(5X,''Parallelizing q-star over'',I2,'' images'')') nimage
  if (nimage > 7) write(stdout,'(5X,''ATTENTION: optimal number of images is 7'')')

  call start_clock('GIPAW')
  call gipaw_readin()
  call check_stop_init()

  io_level = 1
 
  ! read ground state wavefunctions
  cell_factor = 1.1d0
  call read_file
#ifdef __MPI
  use_para_diag = check_para_diag(nbnd)
#else
  use_para_diag = .false.
#endif

  call gipaw_openfil

  if (gamma_only) call errore ('gipaw_main', 'Cannot run GIPAW with gamma_only == .true. ', 1)
  if ((twfcollect .eqv. .false.)  .and. (nproc_pool_file /= nproc_pool)) &
    call errore('gipaw_main', 'Different number of CPU/pool. Set wf_collect=.true. in SCF', 1)
#ifdef __BANDS
  if (nbgrp > 1 .and. (twfcollect .eqv. .false.)) &
    call errore('gipaw_main', 'Cannot use band-parallelization without wf_collect in SCF', 1)
#endif
  if (noncolin) call errore('gipaw_main', 'non-collinear not supported yet', 1)

  call gipaw_allocate()
  call gipaw_setup()
  call gipaw_summary()
  
  ! convert q_gipaw into units of tpiba
  q_gipaw = q_gipaw / tpiba
  
  ! image initialization
  if (nimage >= 1) then
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
  
  ! print timings and stop the code
  CALL gipaw_closefil
  CALL print_clock_gipaw()
  CALL stop_code( .TRUE. )
  
  STOP
  
END PROGRAM gipaw_main

