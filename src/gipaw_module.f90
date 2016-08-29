!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE gipaw_module
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for GIPAW calculations
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : a0_to_cm => bohr_radius_cm
  USE parameters, ONLY : npk, ntypx, lmaxx
  
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER:: nbrx=2*(lmaxx+1)**2  ! max number of beta functions

  ! Physical constants:
  ! fine structure constant alpha (c = 1/alpha)
  REAL(DP), PARAMETER :: alpha = 1.0_dp / 137.03599911_dp
  
  ! avogadro number
  REAL(DP), PARAMETER :: avogadro = 6.022142e23_dp
  
  ! g_e and gprime
  REAL(DP), PARAMETER :: g_e = 2.0023192778_DP
  REAL(DP), PARAMETER :: gprime = 2.d0 * (g_e - 1.d0)
 
  ! rydberg to Hartree
  REAL(DP), PARAMETER :: ry2ha = 0.5_DP
 
  ! number of occupied bands at each k-point
  INTEGER :: nbnd_occ(npk)
  
  ! alpha shift of the projector on the valence wfcs
  REAL(DP) :: alpha_pv
  
  ! eigenvalues and eigenfunctions at k+q
  REAL(DP), ALLOCATABLE :: etq(:,:)
  COMPLEX(DP), ALLOCATABLE :: evq(:,:)

  ! convergence threshold for diagonalizationa and greenfunction
  REAL(DP) :: conv_threshold
  
  ! q for the perturbation (in bohrradius^{-1})
  REAL(DP) :: q_gipaw
  
  ! q for the EFG
  REAL(DP) :: q_efg ( ntypx )

  ! restart mode: 'from_scratch' or 'restart'
  CHARACTER(80) :: restart_mode
  
  ! verbosity
  INTEGER :: iverbosity
 
  ! diagonalization method
  INTEGER :: isolve
 
  ! job: nmr, g_tensor, efg, hyperfine
  CHARACTER(80) :: job
 
  ! core-relax method to calculate change of XC
  INTEGER :: core_relax_method = 1
 
  ! format for a rank-2 tensor
  CHARACTER(*), PARAMETER :: tens_fmt = '(3(5X,3(F14.4,2X)/))'
  
  ! for plotting the induced current and induced field
  CHARACTER(80) :: filcurr, filfield, filnics

  ! macroscopic shape for the NMR
  LOGICAL :: use_nmr_macroscopic_shape
  REAL(DP) :: nmr_macroscopic_shape(3,3)

  ! contribution to NMR chemical shift due to core contribution
  REAL(dp) :: nmr_shift_core(ntypx)
  
  ! parametres for hyper-fine interaction
  CHARACTER(10) :: hfi_input_unit
  CHARACTER(10) :: hfi_output_unit
  REAL(dp) :: hfi_nuclear_g_factor(ntypx)
  LOGICAL :: hfi_via_reconstruction_only = .false. !UWG: speed up for large systems 

  
  REAL(dp) :: rc(ntypx,0:lmaxx)
  COMPLEX(dp), ALLOCATABLE :: paw_becp2(:,:), paw_becp3(:,:)
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: lx, ly, lz
  REAL(dp), ALLOCATABLE :: radial_integral_paramagnetic(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_diamagnetic(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_paramagnetic_so(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_diamagnetic_so(:,:,:)
  REAL(dp), ALLOCATABLE :: radial_integral_rmc(:,:,:)

  LOGICAL :: pawproj( ntypx ) ! EMINE: if true paw projectors will be used instead of GIPAW ones
  REAL(dp) :: r_rand ! EMINE: determines the randimization range used in
                     ! compute_u_kq routine. read from input. change it for debugging only

  INTEGER, PARAMETER :: iumagres = 57

#ifdef __BANDS
  INTEGER :: ibnd_start, ibnd_end
#endif
!-----------------------------------------------------------------------
END MODULE gipaw_module
!-----------------------------------------------------------------------
