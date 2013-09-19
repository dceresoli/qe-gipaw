!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE write_tensor_field(name, ispin, field)
  !-----------------------------------------------------------------------
  !
  ! ... write the tensor field in VTK format
  !
  USE kinds,                       ONLY : dp
  USE io_global,                   ONLY : stdout, ionode
  USE cell_base,                   ONLY : at, alat
  USE ions_base,                   ONLY : nat, tau, atm, ityp
  USE fft_base,                    ONLY : dfftp, grid_gather
  USE pwcom
  USE gipaw_module
  !--------------------------------------------------------------------
  implicit none
  character*(*) name
  integer :: ispin
  real(dp) :: field(dfftp%nnr,3,3)
  !--------------------------------------------------------------------
  integer, parameter :: ounit = 48
  character*80 :: fname
  integer :: ios, ipol
  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: aux(:,:,:)
#ifdef __MPI
  integer :: i, j
#endif

 allocate( aux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3,3) )

#ifdef __MPI
 ! gather the data 
  do i = 1, 3
    do j = 1, 3
      call grid_gather(field(:,i,j), aux(:,i,j))
    enddo
  enddo
#else
  aux = field
#endif

  if (ionode) then
    do ipol = 1, 3
      ! form the name of the output file
      if (nspin == 1) then
        fname = trim(name)//'_'
      elseif (nspin == 2 .and. ispin == 1) then
        fname = trim(name)//'_UP_'
      elseif (ispin == 2 .and. ispin == 2) then
        fname = trim(name)//'_DW_'
      endif

      if (ipol == 1) fname = trim(fname)//'X.vtk'
      if (ipol == 2) fname = trim(fname)//'Y.vtk'
      if (ipol == 3) fname = trim(fname)//'Z.vtk'
      write(stdout, '(5X,''write_tensor_field: '',A40)') fname

      open(unit=ounit, file=fname, iostat=ios, form='formatted', status='unknown')
      if (ios /= 0) call errore('write_tensor_field', 'error opening '//fname, ounit)

      call vtk_vector_3d(aux(:,:,ipol), dfftp%nr1, dfftp%nr2, dfftp%nr3, at, alat, ounit)

      close(unit=ounit)
    enddo
  endif

  deallocate(aux)

end subroutine write_tensor_field


!-------------------------------------------------------------------
! this routine writes a 3D vector field in VTK format
!-------------------------------------------------------------------
subroutine vtk_vector_3d(vin, nr1, nr2, nr3, at, alat, ounit)
  USE kinds,           only : dp
  USE constants,       only : bohr_radius_angs
  USE fft_base,        only : dfftp

  implicit none
  integer, intent(in) ::  nr1, nr2, nr3, ounit
  real(dp), intent(in) :: at(3,3), alat
  real(dp), intent(in) :: vin(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,3)  

  integer :: i1, i2, i3, bla
  real(dp) :: x(3)

  ! header
  write(ounit,'(A)') '# vtk DataFile Version 2.0'
  write(ounit,'(A)') 'created by qe-gipaw'
  write(ounit,'(A)') 'ASCII'
  write(ounit,'(A)') 'DATASET STRUCTURED_GRID'
  write(ounit,'(A,3I5)') 'DIMENSIONS', nr1, nr2, nr3
  write(ounit,'(A,I10,4X,A)') 'POINTS', nr1*nr2*nr3, 'float'

  ! point coordinates
  do i3 = 1, nr3
    do i2 = 1, nr2
      do i1 = 1, nr1
        ! coordinate in angstrom
        x(1) = dble(i1-1)/dble(nr1)     
        x(2) = dble(i2-1)/dble(nr2)     
        x(3) = dble(i3-1)/dble(nr3)
        ! crystal to cartesian
        call cryst_to_cart (1, x, at, 1)
        x = x * alat * BOHR_RADIUS_ANGS
        write(ounit,'(3F15.8)') x
      enddo
    enddo
  enddo

  ! vectors
  write(ounit,'(A,I10)') 'POINT_DATA', nr1*nr2*nr3
  write(ounit,'(A)') 'VECTORS vectors float'
  do i3 = 1, nr3
    do i2 = 1, nr2
      do i1 = 1, nr1
        bla = i1 + (i2-1)*dfftp%nr1x + (i3-1)*dfftp%nr1x*dfftp%nr2x
        write(ounit,'(3E15.8)') vin(bla,1:3)
      enddo
    enddo
  enddo

end subroutine vtk_vector_3d


!-----------------------------------------------------------------------
SUBROUTINE write_nics(filename, field)
  !-----------------------------------------------------------------------
  !
  ! ... write the NICS in PP postproc format
  !
  USE kinds,           ONLY : dp
  USE io_global,       ONLY : stdout, ionode
  USE fft_base,        ONLY : dfftp, grid_gather
  USE gvect,           ONLY : gcutm
  USE cell_base,       ONLY : at, alat, tpiba2, omega, ibrav, celldm
  USE ions_base,       ONLY : zv, ntyp => nsp, nat, ityp, atm, tau
  USE wvfct,           ONLY : ecutwfc
  USE pwcom
  USE gipaw_module
  !--------------------------------------------------------------------
  implicit none
  character*(*) filename
  real(dp) :: field(dfftp%nnr,3,3,nspin)
  !-- local variables ----------------------------------------------------
  character(75), parameter :: title = 'NICS'
  real(dp), allocatable :: nics(:), aux(:)
  integer :: ispin

  allocate(nics(dfftp%nnr))
  nics = 0.d0
  do ispin = 1,nspin
    nics(:) = nics(:) + (field(:,1,1,ispin) + field(:,2,2,ispin) + &
                         field(:,3,3,ispin))/3.d0
  enddo
  nics = nics * 1d6

  ! gather the data 
  allocate(aux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
#ifdef __MPI
  call grid_gather(nics, aux)
#else
  aux = nics
#endif

  if (ionode) then
      write(stdout, '(5X,''writings NICS to: '',A40)') trim(filename)
      call plot_io(filename, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                   dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, &
                   celldm, at, gcutm, dual, ecutwfc, 100, atm, ityp, zv, &
                   tau, aux, 1)
  endif

  deallocate(aux, nics)

end subroutine write_nics

