! generate the actual f90 file with:
! ifort -E xml_routines_FPP.f90 >xml_routines.f90

!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! XML routines non relying on qexsd. For a small code such as
! GIPAW, qexsd is an overkill. Here I use the FoX library directly
! to input and output XML files

! commodity macros
#define CHECK_IONODE if (.not. ionode) return

! XML output
#define _NE(x) call XML_NewElement(xmlf, #x)
#define _ADDS(x) call XML_AddCharacters(xmlf, trim((x)))
#define _ADDV(x) call XML_AddCharacters(xmlf, (x))
#define _EE(x) call XML_EndElement(xmlf, #x)
#define _OUTS(x) call XML_NewElement(xmlf, #x); call XML_AddCharacters(xmlf, trim((x))); call XML_EndElement(xmlf, #x);
#define _OUTV(x) call XML_NewElement(xmlf, #x); call XML_AddCharacters(xmlf, (x)); call XML_EndElement(xmlf, #x);
#define _ATTR(x,y) call XML_AddAttribute(xmlf, #x, (y))

!-----------------------------------------------------------------------
MODULE xml_routines
!-----------------------------------------------------------------------  
  ! ... This module contains the variables routines to I/O the XML files
  USE kinds,                       ONLY : dp
  USE io_global,                   ONLY : ionode
  USE fox_wxml                     ! to output XML
  USE fox_dom                      ! to parse XML file

  IMPLICIT NONE
  SAVE
  
  type(xmlf_t) :: xmlf             ! xml file handle

  character(5),  parameter :: fmt_name = "QEXSD"
  character(8),  parameter :: fmt_version = "21.07.16"

  PUBLIC :: open_xml_output, close_xml_output
  PUBLIC :: xml_output_generalinfo, xml_output_parallelinfo
  PUBLIC :: xml_output_namelist, xml_output_results
  PUBLIC :: xml_output_status
  PUBLIC :: xml_output_timinginfo
  PUBLIC :: xml_output_closed

  CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE open_xml_output(filename)
  !-----------------------------------------------------------------------
    implicit none
    character(len=*), intent(in) :: filename
    integer :: ierr

    CHECK_IONODE
    call XML_OpenFile(filename, xmlf, replace=.true., pretty_print=.true., namespace=.true., iostat=ierr)
    !if ierr

    ! write namespace
    call XML_DeclareNamespace(xmlf, prefix="xsi", nsURI="http://www.w3.org/2001/XMLSchema-instance")
    call XML_DeclareNamespace(xmlf, prefix="qes", nsURI="http://www.quantum-espresso.org/ns/qes/qes-1.0")
    call XML_DeclareNamespace(xmlf, prefix="gpw", nsURI="http://www.quantum-espresso.org/ns/gpw/qes_gipaw_1.0")
    call XML_NewElement(xmlf, name="gpw:gipaw")
    call XML_AddAttribute(xmlf, name="xsi:schemaLocation", &
                          value="http://www.quantum-espresso.org/ns/gpw/qes_gipaw_1.0 "//&
                                "http://www.quantum-espresso.org/ns/gpw/gpw_202201.xsd")
  END SUBROUTINE open_xml_output


  !-----------------------------------------------------------------------
  SUBROUTINE close_xml_output
  !-----------------------------------------------------------------------
    implicit none
   
    CHECK_IONODE
    call XML_Close(xmlf, empty=.false.)

  END SUBROUTINE close_xml_output


  !-----------------------------------------------------------------------
  SUBROUTINE xml_output_generalinfo
  !-----------------------------------------------------------------------
    USE gipaw_version
    implicit none
    character(9) :: cdate, ctime
    character(60) :: timestamp

    call date_and_tim(cdate, ctime)
    timestamp = 'This run was terminated on:  ' // ctime // ' ' // cdate(1:2) // ' '//cdate(3:5) // ' '// cdate (6:9)

    CHECK_IONODE
    _NE(general_info)
      _NE(xml_format)
        _ATTR(NAME, fmt_name)
        _ATTR(VERSION, fmt_version)
        _ADDV(fmt_name//'_'//fmt_version)
      _EE(xml_format)
      _NE(creator)
        _ATTR(NAME, 'GIPAW')
        _ATTR(VERSION, gipaw_git_revision)
        _ADDV('XML file generated by GIPAW')
      _EE(creator)
      _NE(created)
        _ATTR(DATE, cdate)
        _ATTR(TIME, ctime)
        _ADDS(timestamp)
      _EE(created)
      _NE(job)
      _EE(job)
    _EE(general_info)

  END SUBROUTINE xml_output_generalinfo


  !-----------------------------------------------------------------------
  SUBROUTINE xml_output_parallelinfo
  !-----------------------------------------------------------------------
    USE mp_world,         ONLY : nproc
    USE mp_pools,         ONLY : npool
    USE mp_bands,         ONLY : ntask_groups, nbgrp
    USE laxlib_processors_grid, ONLY : nproc_ortho
    implicit none
    integer :: nthreads=1
#if defined(__OMP) 
    integer, external :: omp_get_max
    nthreads = omp_get_max()
#endif      

    CHECK_IONODE
    _NE(parallel_info)
      _NE(nprocs)
        _ADDV(nproc)
      _EE(nprocs)
      _OUTV(nthreads)
      _NE(ntasks)
        _ADDV(ntask_groups)
      _EE(ntasks)
      _OUTV(nbgrp)
      _OUTV(npool)
      _NE(ndiag)
        _ADDV(nproc_ortho)
      _EE(ndiag)
    _EE(parallel_info)

  END SUBROUTINE xml_output_parallelinfo


  !-----------------------------------------------------------------------
  SUBROUTINE xml_output_namelist
  !-----------------------------------------------------------------------
    USE io_files,        ONLY : prefix, tmp_dir
    USE uspp_data,       ONLY : spline_ps
    USE cell_base,       ONLY : tpiba
    USE ions_base,       ONLY : nat
    USE gipaw_module     ! to access internal variables
    implicit none

    CHECK_IONODE
    _NE(input)
      _OUTS(job)
      _OUTS(prefix)
      _OUTS(tmp_dir)
      _OUTV(conv_threshold)
      _OUTS(restart_mode)
      _NE(q_gipaw)
        _ADDV(q_gipaw*tpiba)
      _EE(q_gipaw)
      _OUTS(verbosity)
      _OUTS(filcurr)
      _OUTS(filfield)
      _OUTS(filnics)
      _OUTV(pawproj)
      _OUTV(use_nmr_macroscopic_shape)
      _NE(nmr_macroscopic_shape)
        _ATTR(rank, "2")
        _ATTR(dims, "3 3")
        _ADDV(nmr_macroscopic_shape)
      _EE(nmr_macroscopic_shape)
      _OUTV(spline_ps)
      _OUTS(diagonalization)
      _NE(q_efg)
        _ATTR(size, nat)
        _ADDV(q_efg)
      _EE(q_efg)
      _OUTV(max_seconds)
      _OUTV(r_rand)
      _OUTS(hfi_output_unit)
      _NE(hfi_nuclear_g_factor)
        _ATTR(size, nat)
        _ADDV(hfi_nuclear_g_factor)
      _EE(hfi_nuclear_g_factor)
      _OUTV(core_relax_method)
      _OUTV(hfi_via_reconstruction_only)
    _EE(input)
  END SUBROUTINE xml_output_namelist


  !-----------------------------------------------------------------------
  SUBROUTINE xml_output_susceptibility
  !-----------------------------------------------------------------------
    USE gipaw_results
    implicit none

    CHECK_IONODE

    _NE(susceptibility_low)
      _ATTR(rank, '2')
      _ATTR(dims, '3 3')
      _ADDV(res_suscept1)
    _EE(susceptibility_low)

    _NE(susceptibility_high)
      _ATTR(rank, '2')
      _ATTR(dims, '3 3')
      _ADDV(res_suscept2)
    _EE(susceptibility_high)

  END SUBROUTINE xml_output_susceptibility


  !-----------------------------------------------------------------------
  SUBROUTINE xml_atom_attributes(na, units, rank, dims)
  !-----------------------------------------------------------------------
    USE gipaw_results
    USE ions_base, ONLY: tau, atm, ityp
    implicit none
    integer, intent(in) :: na
    character(*) :: rank, dims, units

    _ATTR(name, trim(atm(ityp(na))))
    _ATTR(tau, tau(:,na))
    _ATTR(index, na)
    _ATTR(rank, rank)
    _ATTR(dims, dims)
    _ATTR(units, units)

  END SUBROUTINE xml_atom_attributes


  !-----------------------------------------------------------------------
  SUBROUTINE xml_output_results
  !-----------------------------------------------------------------------
    USE gipaw_module     ! to access internal variables
    USE gipaw_results
    USE ions_base, ONLY: nat
    implicit none
    integer :: na

    CHECK_IONODE
    _NE(output)

!!    if (job == 'nmr') then
      call xml_output_susceptibility
      _NE(shielding_tensors)
      do na = 1, nat
        _NE(atom)
          call xml_atom_attributes(na, 'ppm', '2', '3 3')
          _ADDV(res_nmr_sigma(:,:,na))
        _EE(atom)
      enddo
      _EE(shielding_tensors)
!!    endif

!!    if (job == 'g_tensor') then
!!      call xml_output_susceptibility
      _NE(delta_g)
        _ATTR(rank, '2')
        _ATTR(dims, '3 3')
        _ADDV(res_epr_deltag)
      _EE(delta_g)
      _NE(delta_g_paratec)
        _ATTR(rank, '2')
        _ATTR(dims, '3 3')
        _ADDV(res_epr_deltag_paratec)
      _EE(delta_g_paratec)
!!    endif

!!    if (job == 'efg') then
      _NE(electric_field_gradients)
        do na = 1, nat
          _NE(atom)
            call xml_atom_attributes(na, 'MHz', '2', '3 3')
            _ADDV(res_efg(:,:,na))
          _EE(atom)
        enddo
        _EE(electric_field_gradients)
!!    endif

!!    if (job == 'hyperfine') then
      _NE(hyperfine_dipolar)
      do na = 1, nat
        _NE(atom)
          call xml_atom_attributes(na, trim(hfi_output_unit), '2', '3 3')
          _ADDV(res_hfi_dip(:,:,na))
        _EE(atom)
      enddo
      _EE(hyperfine_dipolar)
      _NE(hyperfine_fermi_contact)
      do na = 1, nat
        _NE(atom)
          call xml_atom_attributes(na, trim(hfi_output_unit), '1', '1')
          _ADDV(res_hfi_fc(na))
        _EE(atom)
      enddo
      _EE(hyperfine_fermi_contact)
!!    endif

    _EE(output)

  END SUBROUTINE xml_output_results


  !-----------------------------------------------------------------------
  SUBROUTINE xml_output_status
  !-----------------------------------------------------------------------
    implicit none
    integer :: exit_status = 0

    CHECK_IONODE
    _OUTV(exit_status)

  END SUBROUTINE xml_output_status


  !-----------------------------------------------------------------------
  SUBROUTINE xml_output_closed
  !-----------------------------------------------------------------------
    implicit none
    character(9) :: cdate, ctime

    CHECK_IONODE
    call date_and_tim(cdate, ctime)
    _NE(closed)
       _ATTR(DATE, cdate)
       _ATTR(TIME, ctime)
    _EE(closed)

  END SUBROUTINE xml_output_closed


  !-----------------------------------------------------------------------
  SUBROUTINE xml_output_timinginfo
  !-----------------------------------------------------------------------
    USE mytime,  ONLY : f_wall, f_tcpu
    USE gipaw_module
    implicit none
    real(dp), external :: get_clock
    real(dp) :: wall, cpu

    CHECK_IONODE
    _NE(timing_info)
      _NE(total)
        _ATTR(label, 'GIPAW')
        cpu = f_tcpu()
        wall = f_wall()
        !wall = get_clock('GIPAW')
        _OUTV(cpu)
        _OUTV(wall)
      _EE(total)

    _EE(timing_info)

  END SUBROUTINE xml_output_timinginfo


END MODULE xml_routines

