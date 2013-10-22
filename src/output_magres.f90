!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! output info according to magres-v1 format

!-----------------------------------------------------------------------
SUBROUTINE output_magres_begin(what)
  !-----------------------------------------------------------------------
  !
  ! ... Open file and output <atoms> block
  !  
  USE kinds,                  ONLY : dp 
  USE io_files,               ONLY : prefix
  USE io_global,              ONLY : ionode
  USE constants,              ONLY : bohr_radius_angs
  USE cell_base,              ONLY : alat, at
  USE gipaw_module,           ONLY : iumagres
  USE ions_base,              ONLY : nat, atm, ityp, tau
  IMPLICIT NONE
  character(3), intent(in) :: what  ! 'nmr' or 'efg'
  real(dp) :: lattice(3,3), pos(3)
  character(10) :: atom
  integer :: na, ios

  if (.not. ionode) return
  if (what /= 'nmr' .and. what /= 'efg') &
     call errore('output_magres_begin', what // ' not implmented', 1)

  open(unit=iumagres, file=trim(prefix)//'.'//what// '.magres', &
       status='unknown', form='formatted', iostat=ios)
  if (ios /= 0) call errore('output_magres_begin', 'cannot open file', abs(ios))

  write(iumagres,'(''#$magres-abinitio-v1.0'')')
  write(iumagres,*)

  write(iumagres,'(''<atoms>'')')
  write(iumagres,'(''  units lattice Angstrom'')')
  lattice = alat*at * BOHR_RADIUS_ANGS
  write(iumagres,'(''  lattice '',9(F12.6))') lattice(:,1), lattice(:,2), lattice(:,3)
  write(iumagres,*)
  write(iumagres,'(''  units atom Angstrom'')')
  do na = 1, nat
    atom = atm(ityp(na))
    pos = tau(:,na)*alat
    write(iumagres,'(''  atom '',A2,1X,A5,1X,I4,2X,3(F12.6))') atom(1:2), &
          atm(ityp(na)), na, pos*bohr_radius_angs
  enddo
  write(iumagres,'(''</atoms>'')')
  write(iumagres,*)

END SUBROUTINE output_magres_begin


!-----------------------------------------------------------------------
SUBROUTINE output_magres_nmr(chi_bare_pGv, chi_bare_vGv, sigma_tot)
  !-----------------------------------------------------------------------
  !
  ! ... Output <magres> block
  !  
  USE kinds,                  ONLY : dp
  USE io_files,               ONLY : prefix
  USE io_global,              ONLY : ionode
  USE gipaw_module,           ONLY : iumagres, avogadro, a0_to_cm
  USE ions_base,              ONLY : nat, atm, ityp, tau
  USE gipaw_version
  IMPLICIT NONE  
  real(dp), intent(in) :: chi_bare_pGv(3,3), chi_bare_vGv(3,3)
  real(dp), intent(in) :: sigma_tot(3,3,nat)
  real(dp) :: sus(3,3)
  integer :: na

  if (.not. ionode) return

  sus = (chi_bare_pGv + chi_bare_vGv) / 2.d0 * 1e6_dp * a0_to_cm**3.0_dp * avogadro
  write(iumagres,'(''<magres>'')')
  write(iumagres,'(''  units sus 10^-6.cm^3.mol^-1'')')
  write(iumagres,'(''  sus  '',9(F12.4))') sus
  write(iumagres,*)
  write(iumagres,'(''  units ms ppm'')')
  do na = 1, nat
    write(iumagres,'(''  ms '',A5,1X,I4,2X,9(F12.4))') atm(ityp(na)), na, &
          sigma_tot(:,:,na)*1d6
  enddo
  write(iumagres,'(''</magres>'')')
  write(iumagres,*)

END SUBROUTINE output_magres_nmr


!-----------------------------------------------------------------------
SUBROUTINE output_magres_efg(efg_tot)
  !-----------------------------------------------------------------------
  !
  ! ... Output <magres> block
  !  
  USE kinds,                  ONLY : dp
  USE io_files,               ONLY : prefix
  USE io_global,              ONLY : ionode
  USE gipaw_module,           ONLY : iumagres, avogadro, a0_to_cm
  USE ions_base,              ONLY : nat, atm, ityp, tau
  USE gipaw_version
  IMPLICIT NONE  
  real(dp), intent(in) :: efg_tot(3,3,nat)
  integer :: na

  if (.not. ionode) return

  write(iumagres,'(''<magres>'')')
  write(iumagres,'(''  units efg au'')')
  do na = 1, nat
    write(iumagres,'(''  efg '',A5,1X,I4,2X,9(F12.4))') atm(ityp(na)), na, &
          efg_tot(:,:,na)
  enddo
  write(iumagres,'(''</magres>'')')
  write(iumagres,*)

END SUBROUTINE output_magres_efg


!-----------------------------------------------------------------------
SUBROUTINE output_magres_end
  !-----------------------------------------------------------------------
  !
  ! ... Output <calculation> block and close file
  !  
  USE kinds,                  ONLY : dp
  USE io_files,               ONLY : prefix
  USE io_global,              ONLY : ionode
  USE gipaw_module,           ONLY : iumagres
  USE gipaw_version
  IMPLICIT NONE  

  if (.not. ionode) return

  write(iumagres,'(''<calculation>'')')
  write(iumagres,'(''  calc_name          '',A)') trim(prefix)
  write(iumagres,'(''  calc_code          QE-GIPAW 5.x'')')
  write(iumagres,'(''  calc_code_version  svn'',A)') trim(gipaw_svn_revision)
  write(iumagres,'(''</calculation>'')')
  close(iumagres, status='keep')

END SUBROUTINE output_magres_end


