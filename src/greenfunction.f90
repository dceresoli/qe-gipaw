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
  ! ...
  ! ... [H^(0) - alpha P_i - E_i^(0)] |Gpsi_i> = - Q_i H^(1)|psi_i^(0)>
  ! ...
  ! ... We use Hartree atomic units; since G is the inverse of an
  ! ... energy: G => G / ryd_to_hartree
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout  
  USE becmod,                      ONLY : bec_type, becp, calbec, &
                                          allocate_bec_type, deallocate_bec_type
  USE wavefunctions_module,        ONLY : evc
  USE noncollin_module,            ONLY : npol
  USE pwcom,                       ONLY : ef
  USE wvfct,                       ONLY : nbnd, et, npw, npwx, igk, g2kin
  USE gvect,                       ONLY : g
  USE uspp,                        ONLY : nkb, vkb
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE mp_pools,                    ONLY : intra_pool_comm, inter_pool_comm
  USE mp_bands,                    ONLY : inter_bgrp_comm
  USE mp,                          ONLY : mp_sum
  USE ldaU,                        ONLY : lda_plus_u, wfcU
  USE io_files,                    ONLY : iunhub, nwordwfcU
  USE buffers,                     ONLY : get_buffer
  USE cell_base,                   ONLY : tpiba
  USE klist,                       ONLY : lgauss, xk, degauss, ngauss
  USE gipaw_module
#ifdef __BANDS
  USE mp_bands,                    ONLY : intra_bgrp_comm
#endif

  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nbnd)  ! psi is H1*psi and is changed on output!!!
  COMPLEX(DP), INTENT(OUT) :: g_psi(npwx,nbnd)
  REAL(DP) :: q(3)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  complex(dp), allocatable :: ps(:,:), work (:)
  real(dp), allocatable :: h_diag (:,:), eprec (:)
  real(dp) :: anorm, thresh, gk(3), dxk(3)
  integer :: ibnd, jbnd, ig, lter
  logical :: conv_root, q_is_zero
  complex(dp), external :: zdotc
  real(dp), external :: wgauss, w0gauss
  real(dp) :: wg1, w0g, wgp, wwg, deltae, theta
  external ch_psi_all, cg_psi
 
  ! start clock
  call start_clock ('greenf')

  ! allocate memory
  allocate (work(npwx), ps(nbnd,nbnd), h_diag(npwx,nbnd), eprec(nbnd))
  call allocate_bec_type(nkb, nbnd, becp)

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

  if (lgauss) then
     ! metallic case
#ifdef __BANDS
     CALL zgemm('C', 'N', nbnd, ibnd_end-ibnd_start+1, npw, &
                 (1.d0,0.d0), evq(1,1), npwx, psi(1,ibnd_start), npwx, (0.d0,0.d0), &
                 ps(1,ibnd_start), nbnd)
#else
     CALL zgemm('C', 'N', nbnd, nbnd_occ (ik), npw, &
                (1.d0,0.d0), evq(1,1), npwx, psi(1,1), npwx, (0.d0,0.d0), &
                ps(1,1), nbnd)
#endif   
     do ibnd = 1, nbnd_occ(ik)
        wg1 = wgauss ((ef-et(ibnd,ik)) / degauss, ngauss)
        w0g = w0gauss((ef-et(ibnd,ik)) / degauss, ngauss) / degauss
        do jbnd = 1, nbnd
           wgp = wgauss ( (ef - etq(jbnd,ik)) / degauss, ngauss)
           deltae = etq(jbnd,ik) - et(ibnd,ik)
           theta = wgauss (deltae / degauss, 0)
           wwg = wg1 * (1.d0 - theta) + wgp * theta
           if (jbnd <= nbnd_occ(ik)) then
              if (abs (deltae) > 1d-5) then
                 wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
              else
                 ! if the two energies are too close takes the limit of the 0/0 ratio
                 wwg = wwg - alpha_pv * theta * w0g
              endif
           endif
           ps(jbnd,ibnd) = wwg * ps(jbnd,ibnd)
        enddo
        call dscal(2*npw, wg1, psi(1,ibnd), 1)
     enddo

  else
     ! insulators
#ifdef __BANDS
     CALL zgemm('C', 'N', nbnd_occ (ik), ibnd_end-ibnd_start+1, npw, &
                (1.d0,0.d0), evq(1,1), npwx, psi(1,ibnd_start), npwx, (0.d0,0.d0), &
                ps(1,ibnd_start), nbnd)
#else
     CALL zgemm('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
                (1.d0,0.d0), evq(1,1), npwx, psi(1,1), npwx, (0.d0,0.d0), &
                ps(1,1), nbnd)
#endif   
  endif

#ifdef __MPI
#ifdef __BANDS
  call mp_sum(ps, intra_bgrp_comm)
#else
  call mp_sum(ps, intra_pool_comm)
#endif
#endif


  ! this is the case with overlap (ultrasoft)
  ! g_psi is used as work space to store S|evq>
  ! |psi> = -(|psi> - S|evq><evq|psi>)
#ifdef __BANDS
  CALL calbec_bands (npwx, npw, nkb, vkb, evq, becp%k, nbnd_occ(ik), ibnd_start, ibnd_end)
  CALL s_psi_bands (npwx, npw, nbnd_occ(ik), evq, g_psi, ibnd_start, ibnd_end)
#else
  CALL calbec (npw, vkb, evq, becp)
  if (lgauss) then 
     CALL s_psi (npwx, npw, nbnd, evq, g_psi)
  else
     CALL s_psi (npwx, npw, nbnd_occ(ik), evq, g_psi)
  endif
#endif

#if defined(__MPI) && defined(__BANDS)
  ! replicate wfc
  call mp_sum(g_psi, inter_bgrp_comm)
#endif

  if (lgauss) then
     ! metallic case
#ifdef __BANDS
     CALL zgemm( 'N', 'N', npw, ibnd_end-ibnd_start+1, nbnd, &
          (1.d0,0.d0), g_psi(1,1), npwx, ps(1,ibnd_start), nbnd, (-1.d0,0.d0), &
          psi(1,ibnd_start), npwx )
#else
     CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd, &
          (1.d0,0.d0), g_psi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
          psi(1,1), npwx )
#endif

  else
     ! insulators
#ifdef __BANDS
     CALL zgemm( 'N', 'N', npw, ibnd_end-ibnd_start+1, nbnd_occ(ik), &
          (1.d0,0.d0), g_psi(1,1), npwx, ps(1,ibnd_start), nbnd, (-1.d0,0.d0), &
          psi(1,ibnd_start), npwx )
#else
     CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
          (1.d0,0.d0), g_psi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
          psi(1,1), npwx )
#endif
endif

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
  do ibnd = 1, nbnd
#endif
     do ig = 1, npw
        work (ig) = g2kin (ig) * evq (ig, ibnd)
     enddo
     eprec (ibnd) = 1.35d0 * zdotc (npw, evq (1, ibnd), 1, work, 1)
  enddo
#ifdef __MPI
#ifdef __BANDS
  call mp_sum ( eprec, intra_bgrp_comm )
#else
  call mp_sum ( eprec, intra_pool_comm )
#endif
#endif
  h_diag = 0.d0
#ifdef __BANDS
  do ibnd = ibnd_start, ibnd_end
#else
  do ibnd = 1, nbnd
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

#ifdef __BANDS
  call calbec_bands (npwx, npw, nkb, vkb, psi, becp%k, nbnd, ibnd_start, ibnd_end)
#else
  call calbec (npw, vkb, psi, becp, nbnd)
#endif

  if (lda_plus_u) then
    if (q_is_zero) then
      call get_buffer(wfcU, nwordwfcU, iunhub, ik)
    else
      call orthoatwfc1(ik)
    endif
  endif

  ! initial guess
  g_psi(:,:) = (0.d0, 0.d0)

  ! solve linear system  
  conv_root = .true.
  call cgsolve_all (ch_psi_all, cg_psi, et(1,ik), psi, g_psi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
       nbnd_occ(ik), npol )

#if defined(__MPI) && defined(__BANDS)
  ! replicate wfc
  call mp_sum(g_psi, inter_bgrp_comm)
#endif

  if (iverbosity > 20) &
    write(stdout, '(5X,''cgsolve_all iterations '',I3,4X,''anorm='',E12.2)')  lter, anorm

  ! convert to Hartree
  g_psi(:,:) = g_psi(:,:) / ryd_to_hartree

  call flush_unit( stdout )
  call stop_clock('greenf')
 
  ! free memory
  deallocate (work, h_diag, eprec, ps)
  call deallocate_bec_type (becp)

END SUBROUTINE greenfunction
