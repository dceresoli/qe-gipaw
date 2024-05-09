! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! The velocity operator is composed of three terms:
! v(\epsilon) = (1/i)[r, H - \epsilon S] = p + (1/i)[r, V_NL] - \epsilon (1/i)[r, S]
!
! the three terms are calculated in the following routines 
! apply_p      => apply p to the wave functions
! apply_vel_NL => apply (1/i)[r,V_NL] or (1/i)[r,S] to the wave functions
!
! Finally, the apply_vel subroutine is a driver that applies the velocity operator
!
!-----------------------------------------------------------------------
SUBROUTINE apply_p(psi, p_psi, ik, ipol, q)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the kinetic part of the velocity operator
  ! ... |p_psi> = (G+k+q/2)_{ipol} |psi>
  !  
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, igk_k, ngk
  USE wvfct,                ONLY : nbnd, npwx
  USE gipaw_module,         ONLY : nbnd_occ
  USE gvect,                ONLY : g
  USE cell_base,            ONLY : tpiba
#ifdef __BANDS
  USE mp_bands,             ONLY : inter_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE gipaw_module,         ONLY : ibnd_start, ibnd_end
#endif

  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik               ! k-point
  INTEGER, INTENT(IN) :: ipol             ! cartesian direction (1..3)
  REAL(DP), INTENT(IN) :: q(3)
  COMPLEX(DP), INTENT(IN) :: psi(npwx,nbnd)
  COMPLEX(DP), INTENT(OUT) :: p_psi(npwx,nbnd)

  !-- local variables ----------------------------------------------------
  REAL(DP) :: gk
  INTEGER :: ig, ibnd
  INTEGER :: npw

  do ibnd = 1, nbnd_occ(ik)
    npw = ngk(ik)
    do ig = 1, npw
      gk = xk(ipol,ik) + g(ipol,igk_k(ig,ik)) + q(ipol)
      p_psi(ig,ibnd) = p_psi(ig,ibnd) + gk * tpiba * psi(ig,ibnd)
    enddo
  enddo

END SUBROUTINE apply_p

 
!-----------------------------------------------------------------------
SUBROUTINE apply_vel_NL(what, psi, vel_psi, ik, ipol, q)
  !-----------------------------------------------------------------------
  !
  ! ... Apply (what = 'V') the non-local part of the velocity operator:
  ! ...   (1/i)[r,V_NL] => dV^{NL}_{k+q,k}/dk
  ! ... here we use Hartree atomic units, so that:
  ! ...   V^{NL} => V^{NL} * ryd_to_hartree
  ! ...
  ! ... or (what = 'S') the ultrasoft contribution:
  ! ...   (1/i)[r,S] => dS_{k+q,k}/dk
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, igk_k, ngk
  USE wvfct,                ONLY : nbnd, npwx, current_k
  USE becmod,               ONLY : bec_type, becp, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  USE uspp,                 ONLY : nkb, vkb, qq_at
  USE cell_base,            ONLY : tpiba
  USE gipaw_module,         ONLY : q_gipaw, nbnd_occ
  USE ldaU,                 ONLY : lda_plus_U
  USE lsda_mod,             ONLY : current_spin, lsda, isk, nspin
#ifdef __BANDS
  USE mp_bands,             ONLY : inter_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE gipaw_module,         ONLY : ibnd_start, ibnd_end
#endif
  !-- paramters ----------------------------------------------------------
  IMPLICIT NONE
  CHARACTER, INTENT(IN) :: what     ! 'S' or 'V'
  INTEGER, INTENT(IN) :: ipol       ! cartesian direction (1..3)
  INTEGER, INTENT(IN) :: ik         ! k-point
  COMPLEX(DP), INTENT(IN) :: psi(npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: vel_psi(npwx,nbnd)
  REAL(DP), INTENT(IN) :: q(3)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  complex(dp), allocatable :: aux(:,:), vkb_save(:,:)
#ifdef __BANDS
  complex(dp), allocatable :: aux2(:,:)
#endif
  real(dp) :: dk, dxk(3)
  integer :: isign
  logical :: q_is_zero
  integer :: npw

  ! if no projectors, return
  if (nkb == 0) return

  ! set dk
  dk = q_gipaw/2.d0

  ! check if |q| is zero
  q_is_zero = .false.
  if (sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3)) < 1d-8) q_is_zero = .true.

  ! initialization
  npw = ngk(ik)
  current_k = ik
  if (lsda) current_spin = isk(ik)

  ! allocate temporary arrays, save old NL-potential
  allocate(aux(npwx,nbnd), vkb_save(npwx,nkb))
  vkb_save = vkb

#ifdef __BANDS
  allocate(aux2(npwx,nbnd))
  aux2 = (0.0d0, 0.0d0)
#endif

  !====================================================================
  ! compute (1/2|dk|) ( V^{NL}_{k+dk+q,k+dk} |psi> - 
  !                     V^{NL}_{k-dk+q,k-dk} |psi> )
  ! or the same, with S.
  !====================================================================
  do isign = -1,1,2
      dxk(:) = xk(:,ik)
      dxk(ipol) = dxk(ipol) + isign * dk     ! k \pm dk

      ! compute <\beta(k \pm dk)| and project on |psi>
      call init_us_2_no_phase(npw, igk_k(1,ik), dxk, vkb, .false.)  ! TODO: .true. if on GPU
      call allocate_bec_type(nkb, nbnd_occ(ik), becp)
#ifdef __BANDS
      call calbec_bands (npwx, npw, nkb, vkb, psi, becp%k, nbnd_occ(ik), ibnd_start, ibnd_end)
#else
      call calbec (npw, vkb, psi, becp, nbnd_occ(ik))
#endif

      ! |q|!=0 => compute |\beta(k \pm dk + q)>
      if (.not. q_is_zero) then
          dxk(:) = dxk(:) + q(:)
          call init_us_2_no_phase(npw, igk_k(1,ik), dxk, vkb, .false.) ! TODO: .true. if on GPU
      endif

      aux = (0.d0,0.d0)
      if (what == 'V' .or. what == 'v') then
          ! apply |\beta(k \pm dk+q)>D<\beta(k \pm dk)| to |psi>
#ifdef __BANDS
          call add_vuspsi_bands(npwx, npw, nbnd_occ(ik), aux, ibnd_start, ibnd_end)
          !! specialized Hubbard term missing
          aux2(:,ibnd_start:ibnd_end) = aux2(:,ibnd_start:ibnd_end) + dble(isign) * ryd_to_hartree * &
                                        aux(:,ibnd_start:ibnd_end)/(2.d0*dk*tpiba)
#else
          call add_vuspsi(npwx, npw, nbnd_occ(ik), aux)
          call deallocate_bec_type(becp)

          if (lda_plus_U) then
             call orthoatwfc1(ik)
             call vhpsi(npwx, npw, nbnd_occ(ik), psi, aux)
          endif
          vel_psi = vel_psi + dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)
#endif

      elseif (what == 'S' .or. what == 's') then
          ! apply |\beta(k \pm dk+q)>S<\beta(k \pm dk)| to |psi>
#ifdef __BANDS
          call s_psi_bands(npwx, npw, nbnd_occ(ik), psi, aux, ibnd_start, ibnd_end)
          aux2(:,ibnd_start:ibnd_end) = aux2(:,ibnd_start:ibnd_end) + dble(isign) * aux(:,ibnd_start:ibnd_end)/(2.d0*dk*tpiba)
#else
          call s_psi(npwx, npw, nbnd_occ(ik), psi, aux)
          call deallocate_bec_type(becp)
          vel_psi = vel_psi + dble(isign) * aux/(2.d0*dk*tpiba)
#endif
      else
          call errore('apply_vel_NL', '''what'' parameter has the wrong value', 1)
      endif 
  enddo

#if defined(__MPI) && defined(__BANDS)  
  call mp_sum(aux2, inter_bgrp_comm)
  vel_psi = vel_psi + aux2
#endif

  ! restore NL-potential at k
  vkb = vkb_save
  
  ! free memory
  deallocate(aux, vkb_save)
#ifdef __BANDS
  deallocate(aux2)
#endif

END SUBROUTINE apply_vel_NL


!-----------------------------------------------------------------------
SUBROUTINE apply_vel(psi, vel_psi, ik, ipol, q)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the velocity operator
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : nbnd, npwx, et 
  USE uspp,                 ONLY : okvan
  USE gipaw_module,         ONLY : nbnd_occ
#ifdef __BANDS
  USE gipaw_module,         ONLY : ibnd_start, ibnd_end
#endif

  !-- paramters ----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ipol       ! cartesian direction (1..3)
  INTEGER, INTENT(IN) :: ik         ! k-point
  COMPLEX(DP), INTENT(IN) :: psi(npwx,nbnd)
  COMPLEX(DP), INTENT(OUT) :: vel_psi(npwx,nbnd)
  REAL(DP), INTENT(IN) :: q(3)

  !-- local variables ----------------------------------------------------
  integer :: ibnd
  real(dp), parameter :: ryd_to_hartree = 0.5d0

  call start_clock('apply_vel')
  vel_psi = (0.d0,0.d0)

  if (okvan) then
      call apply_vel_NL('S', psi, vel_psi, ik, ipol, q)
#ifdef __BANDS
      do ibnd = ibnd_start, ibnd_end
#else
      do ibnd = 1, nbnd_occ(ik)
#endif
          vel_psi(1:npwx,ibnd) = -et(ibnd,ik) * ryd_to_hartree * vel_psi(1:npwx,ibnd)
      enddo
  endif

  call apply_vel_NL('V', psi, vel_psi, ik, ipol, q)

  call apply_p(psi, vel_psi, ik, ipol, q)

  call stop_clock('apply_vel')

END SUBROUTINE apply_vel

