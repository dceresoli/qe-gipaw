! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE greenfunction(ik, psi, g_psi, q)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the Green function operator
  ! ... (formula here)
  ! ... We use Hartree atomic units; since G is the inverse of an
  ! ... energy: G => G / ryd_to_hartree
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout  
  USE becmod,                      ONLY : bec_type, becp, calbec, &
                                          allocate_bec_type, deallocate_bec_type, calbec_new
  USE wavefunctions_module,        ONLY : evc
  USE noncollin_module,            ONLY : npol
  USE pwcom
  USE uspp,                        ONLY : nkb, vkb
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE gipaw_module
  USE mp_global,                   ONLY : intra_bgrp_comm, inter_bgrp_comm,intra_pool_comm, inter_pool_comm
  USE mp,                          ONLY : mp_sum
  USE mp_image_global_module,      ONLY : inter_image_comm, nimage


  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nbnd)  ! psi is changed on output!!!
  COMPLEX(DP), INTENT(OUT) :: g_psi(npwx,nbnd)
  REAL(DP) :: q(3)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  complex(dp), allocatable :: ps(:,:), work (:)
  real(dp), allocatable :: h_diag (:,:), eprec (:)
  real(dp) :: anorm, thresh, gk(3), dxk(3)
  integer :: ibnd, ig, lter
  logical :: conv_root, q_is_zero
  complex(dp), external :: zdotc
  external ch_psi_all, cg_psi
 
  ! start clock
  call start_clock ('greenf')

  ! allocate memory
  allocate (work(npwx), ps(nbnd,nbnd), h_diag(npwx,nbnd), eprec(nbnd))
  call allocate_bec_type ( nkb, nbnd, becp)

  ! check if |q| is zero
  q_is_zero = .false.
  if (sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3)) < 1d-8) q_is_zero = .true.  

  ! evq is already calculated in compute_u_kq.f90
 if (q_is_zero) evq(:,:) = evc(:,:)

  !====================================================================
  ! apply -Q_{k+q} to the r.h.s.
  !====================================================================
  ! project on <evq|: ps(i,j) = <evq(i)|psi(j)>
  ps = (0.d0,0.d0)
#ifdef __BANDS
  CALL zgemm('C', 'N', nbnd_occ (ik), ibnd_end-ibnd_start+1, npw, &
             (1.d0,0.d0), evq(1,1), npwx, psi(1,ibnd_start), npwx, (0.d0,0.d0), &
             ps(1,ibnd_start), nbnd)
#else
  CALL zgemm('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
             (1.d0,0.d0), evq(1,1), npwx, psi(1,1), npwx, (0.d0,0.d0), &
             ps(1,1), nbnd)
#endif   
#ifdef __PARA
#ifdef __BANDS
  call mp_sum(ps,intra_bgrp_comm)
#else
  call mp_sum(ps,intra_pool_comm)
#endif
#endif

  ! this is the case with overlap (ultrasoft)
  ! g_psi is used as work space to store S|evq>
  ! |psi> = -(|psi> - S|evq><evq|psi>)
#ifdef __BANDS
  CALL calbec_new (npw, vkb, evq, becp, nbnd_occ(ik), ibnd_start, ibnd_end)
  CALL s_psi_new (npwx, npw, nbnd_occ(ik), evq, g_psi, ibnd_start, ibnd_end)
#else
  CALL calbec (npw, vkb, evq, becp)
  CALL s_psi (npwx, npw, nbnd_occ(ik), evq, g_psi)
#endif
#ifdef __PARA 
#ifdef __BANDS
  call mp_sum(g_psi, inter_bgrp_comm)
#else
  call mp_sum(g_psi, inter_bgrp_comm)
#endif
#endif
#ifdef __BANDS
  CALL zgemm( 'N', 'N', npw, ibnd_end-ibnd_start+1, nbnd_occ(ik), &
       (1.d0,0.d0), g_psi(1,1), npwx, ps(1,ibnd_start), nbnd, (-1.d0,0.d0), &
       psi(1,ibnd_start), npwx )
#else
  CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
       (1.d0,0.d0), g_psi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
       psi(1,1), npwx )
#endif
  !! this is the old code for norm-conserving:
  !! |psi> = -(1 - |evq><evq|) |psi>
  !!CALL zgemm('N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
  !!           (1.d0,0.d0), evq(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
  !!           psi(1,1), npwx)


  !====================================================================
  ! solve the linear system (apply G_{k+q})
  !====================================================================
  ! convergence treshold
  thresh = sqrt(conv_threshold)   ! sqrt(of that of PARATEC)

  ! use the hamiltonian at k+q
  do ig = 1, npw
    gk(1) = (xk(1,ik) + g(1,igk(ig)) + q(1)) * tpiba
    gk(2) = (xk(2,ik) + g(2,igk(ig)) + q(2)) * tpiba
    gk(3) = (xk(3,ik) + g(3,igk(ig)) + q(3)) * tpiba
    g2kin (ig) = gk(1)**2 + gk(2)**2 + gk(3)**2
  enddo

  ! preconditioning of the linear system
  work = (0.d0,0.d0)
#ifdef __BANDS
  do ibnd = ibnd_start, ibnd_end
#else
  do ibnd = 1, nbnd_occ (ik)
#endif
     do ig = 1, npw
        work (ig) = g2kin (ig) * evq (ig, ibnd)
     enddo
     eprec (ibnd) = 1.35d0 * zdotc (npw, evq (1, ibnd), 1, work, 1)
  enddo
#ifdef __PARA
#ifdef __BANDS
  call mp_sum ( eprec( 1:nbnd_occ(ik) ), intra_bgrp_comm )
#else
  call mp_sum ( eprec( 1:nbnd_occ(ik) ), intra_pool_comm )
#endif
#endif
  h_diag = 0.d0
#ifdef __BANDS
  do ibnd = ibnd_start, ibnd_end
#else
  do ibnd = 1, nbnd_occ (ik)
#endif
     do ig = 1, npw
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
     enddo
  enddo

  if (.not. q_is_zero) then
    dxk = xk(:,ik) + q
    call init_us_2(npw, igk, dxk, vkb)
  else
    call init_us_2(npw, igk, xk(1,ik), vkb)
  endif
  !it was: call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, psi)
#ifdef __BANDS
  call calbec_new (npw, vkb, psi, becp, nbnd, ibnd_start, ibnd_end)
#else
  call calbec (npw, vkb, psi, becp, nbnd)
#endif    
  ! initial guess
  g_psi(:,:) = (0.d0, 0.d0)
  ! solve linear system  
  conv_root = .true.
  call cgsolve_all (ch_psi_all, cg_psi, et(1,ik), psi, g_psi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
       nbnd_occ(ik), npol )
#ifdef __PARA
#ifdef __BANDS
  call mp_sum(g_psi,inter_bgrp_comm)
#endif
#endif
  !! debug  
  !!write(stdout, '(5X,''cgsolve_all converged in '',I3,'' iterations'')') &
  !!      lter

  if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
       & " linter: root not converged ",e10.3)') &
       ik, ibnd, anorm

  ! convert to Hartree
  g_psi(:,:) = g_psi(:,:) / ryd_to_hartree

  call flush_unit( stdout )
  call stop_clock('greenf')
 
  ! free memory
  deallocate (work, h_diag, eprec, ps)
  call deallocate_bec_type (becp)

END SUBROUTINE greenfunction
