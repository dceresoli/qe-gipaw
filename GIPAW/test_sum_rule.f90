!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE test_f_sum_rule
  !-----------------------------------------------------------------------
  !
  ! ... The generalized f-sum rule:
  ! ... -N_{el} \delta{\alpha,\beta} = "p-G-v" term + "occ" term
  ! ...
  ! ... The "p-G-v" term is:
  ! ... 2 sum_k w_k sum_n^{occ} <u_{nk}| p_{k,\alpha} G_k v_{k,\beta} | u_{nk}>
  ! ... 
  ! ... where: p_k = -i\nabla + k
  ! ...        v_k = -i [r, H_k]  
  ! ...
  ! ... The "occ" term is:
  ! ... -2 sum_k w_k sum_{n,m}^{occ} <u_{nk}|p_\alpha|u_{mk}><u_{mk}|(1/i)[r,S]_\beta|u_{nk}>
  !

  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE cell_base,                   ONLY : at, bg, omega, tpiba, tpiba2
  USE wavefunctions_module,        ONLY : evc
  USE klist,                       ONLY : nks, nkstot, wk, xk, nelec
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                          current_k
  USE lsda_mod,                    ONLY : current_spin, lsda, isk
  USE buffers,                     ONLY : get_buffer
  USE gvect,                       ONLY : ngm, g, ecutwfc
  USE uspp,                        ONLY : vkb, okvan
  USE gipaw_module,                ONLY : nbnd_occ, iverbosity
  USE mp_global,                   ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                          ONLY : mp_sum
  USE symme,                       ONLY : symmatrix


  !-- local variables ----------------------------------------------------
  IMPLICIT NONE
  complex(dp), allocatable, dimension(:,:,:) :: p_evc, vel_evc, g_vel_evc
  integer :: ik, ipol, jpol, ibnd, jbnd
  real(dp) :: q(3)
  real(dp) :: f_sum1(3,3), f_sum_k1(3,3)
  real(dp) :: f_sum2(3,3), f_sum_k2(3,3)
  complex(dp), external :: zdotc

  call start_clock("f-sum")

  ! allocate memory
  allocate ( p_evc(npwx,nbnd,3),   &
             vel_evc(npwx,nbnd,3), &
             g_vel_evc(npwx,nbnd,3) )

  ! zero the f-sum
  f_sum1(:,:) = 0.d0
  f_sum2(:,:) = 0.d0
  q(:) = 0.d0

  write(stdout, '(5X,''Computing the (generalized) f-sum rule'')')

  !====================================================================
  ! loop over k-points
  !====================================================================
  do ik = 1, nks
    current_k = ik
    current_spin = isk(ik)

    ! initialize at k-point k 
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(:) = g2kin(:) * tpiba2
    call init_us_2(npw,igk,xk(1,ik),vkb)

    ! read wfcs from file
    call get_buffer (evc, nwordwfc, iunwfc, ik)

    q(:) = 0.0_dp

    ! compute p_k|evc>, v_k|evc> and G_k v_k|evc>
    p_evc(:,:,:) = (0.d0,0.d0)
    vel_evc(:,:,:) = (0.d0,0.d0)
    do ipol = 1, 3
      call apply_p(evc, p_evc(1,1,ipol), ik, ipol, q)
      call apply_vel(evc, vel_evc(1,1,ipol), ik, ipol, q)
      call greenfunction(ik, vel_evc(1,1,ipol), g_vel_evc(1,1,ipol), q)
    enddo

    if (okvan) then
        ! reusing array vel_evc to store dS/dk|evc>
        vel_evc(:,:,:) = (0.d0,0.d0)
        do ipol = 1, 3
            call apply_vel_NL('S', evc, vel_evc(1,1,ipol), ik, ipol, q)
        enddo
    endif

    ! k-point contribution to the f-sum rule
    f_sum_k1 = 0.0d0
    f_sum_k2 = 0.0d0

    ! loop over cartesian directions
    do jpol = 1, 3
      do ipol = 1, 3

        do ibnd = 1, nbnd_occ (ik)
          ! this is the "p-G-v" term
          f_sum_k1(ipol,jpol) = f_sum_k1(ipol,jpol) + wg(ibnd,ik) * &
            2.d0 * real(zdotc(npw, p_evc(1,ibnd,ipol), 1, &
                                   g_vel_evc(1,ibnd,jpol), 1))
          ! this is the "occ" term
          if (okvan) then
              do jbnd = 1, nbnd_occ (ik)
                 f_sum_k2(ipol,jpol) = f_sum_k2(ipol,jpol) + wg(ibnd,ik) * &
                   (-1.d0) * zdotc(npw, evc(1,ibnd), 1, p_evc(1,jbnd,ipol), 1) * &
                             zdotc(npw, evc(1,jbnd), 1, vel_evc(1,ibnd,jpol), 1)
              enddo  ! jbnd
          endif
        enddo   ! ibnd

      enddo   ! ipol
    enddo   ! jpol

    if (iverbosity > 10) then
      write(stdout, '(5X,''f-sum rule (ik='',I5,''):'')') ik
      write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum_k1
      write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum_k2
    endif

    f_sum1(:,:) = f_sum1(:,:) + f_sum_k1(:,:)
    f_sum2(:,:) = f_sum2(:,:) + f_sum_k2(:,:)
  enddo   ! ik
#ifdef __PARA
  call mp_sum( f_sum1, intra_pool_comm )
  call mp_sum( f_sum1, inter_pool_comm )
  call mp_sum( f_sum2, intra_pool_comm )
  call mp_sum( f_sum2, inter_pool_comm )
#endif

  write(stdout, '(5X,''f-sum rule (pGV term):'')')
  write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum1
  write(stdout, '(5X,''f-sum rule (occ term):'')')
  write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum2

  call symmatrix(f_sum1)
  call symmatrix(f_sum2)
  write(stdout, '(5X,''f-sum rule (total and symmetrized):'')')
  write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum1 + f_sum2

  deallocate(p_evc, vel_evc, g_vel_evc)

  call stop_clock("f-sum")

END SUBROUTINE test_f_sum_rule


