!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE suscept_crystal
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the f-sum rule, the magnetic susceptibility and
  ! the induced current, using the "crystal method".
  ! 
  ! References:
  !   [1] Phys. Rev. B 63, 245101 (2001)   (norm-conserving GIPAW)
  !   [2] Phys. Rev. B 76, 024401 (2007)   (ultrasoft)
  !   [3] Phys. Rev. B 76, 165122 (2007)   (metallic systems)
  !
  ! Contributors:
  !   D. Ceresoli                        bare susceptibility and current
  !   A. P. Seitsonen and U. Gerstmann   GIPAW contributions
  !   E. Kucukbenli                      Ultrasoft and PAW
  !
  USE kinds,                  ONLY : dp
  USE io_global,              ONLY : stdout
  USE io_files,               ONLY : nwordwfc, iunwfc
  USE cell_base,              ONLY : at, bg, omega, tpiba, tpiba2
  USE wavefunctions_module,   ONLY : evc
  USE klist,                  ONLY : nks, nkstot, wk, xk, nelec
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k
  USE lsda_mod,               ONLY : current_spin, lsda, isk
  USE becmod,                 ONLY : becp, calbec  
  USE symme,                  ONLY : symmatrix
  USE parameters,             ONLY : lmaxx
  USE constants,              ONLY : pi
  USE gvect,                  ONLY : ngm, g, ecutwfc
  USE gsmooth,                ONLY : nrxxs
  USE uspp,                   ONLY : vkb, okvan
  USE lsda_mod,               ONLY : nspin
  USE gipaw_module,           ONLY : j_bare, b_ind, b_ind_r, tens_fmt, &
                                     q_gipaw, iverbosity, alpha, evq, &
                                     avogadro, filcurr, filfield, &
                                     nbnd_occ, a0_to_cm, isolve, &
                                     conv_threshold, job
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp, paw_nkb, paw_recon
  USE ions_base,              ONLY : nat
  USE buffers,                ONLY : get_buffer
  USE mp_global,              ONLY : my_pool_id, me_pool, root_pool, &
                                     inter_pool_comm, intra_pool_comm
  USE mp,                     ONLY : mp_sum
  
  !-- local variables ----------------------------------------------------
  IMPLICIT NONE

  ! the following three quantities are for norm-conserving PPs
  complex(dp), allocatable, dimension(:,:,:) :: p_evc         ! p_k|evc>
  complex(dp), allocatable, dimension(:,:,:) :: vel_evc       ! v_{k+q,k}|evc>
  complex(dp), allocatable, dimension(:,:,:) :: G_vel_evc     ! G_{k+q} v_{k+q,k}|evc>
  ! in addition, the following two quantities are for ultrasoft PPs
  complex(dp), allocatable, dimension(:,:,:) :: svel_evc      ! s_{k+q,k}|evc>
  complex(dp), allocatable, dimension(:,:,:) :: u_svel_evc    ! sum|evq><evq|s_{k+q,k}|evc>
  ! temporary working array, same size as evc/evq
  complex(dp), allocatable :: aux(:,:)

  ! f-sum rule: Eq.(A7) of [1]
  real(dp) :: f_sum(3,3)                             ! Eq.(C7) of [2]

  ! Susceptibility (pGv => HH in Paratec, vGv => VV in Paratec)
  real(dp) :: q_pGv(3,3,-1:1), q_vGv(3,3,-1:1)       ! Eq.(65) of [1]
  real(dp) :: f_pGv(3,3,-1:1), f_vGv(3,3,-1:1)       ! Eq.(64) of [1]
  real(dp) :: chi_bare_pGv(3,3), chi_bare_vGv(3,3)   ! Eq.(64) of [1]

  ! GIPAW terms
  real(dp) :: diamagnetic_corr_tensor(3,3,nat)       ! Eq.(58) of [1]/Eq.(43) of [2]
  real(dp) :: paramagnetic_corr_tensor(3,3,nat)      ! Eq.(44) of [1]
  real(dp) :: paramagnetic_corr_tensor_us(3,3,nat)   ! Eq.(41) of [2] ("occ-occ" term)
  real(dp) :: paramagnetic_corr_tensor_aug(3,3,nat)  ! Eq.(30) of [2] (L_R Q_R term)

  ! Contributions to chemical shift
  real(dp) :: sigma_bare(3,3,nat)
  real(dp) :: sigma_diamagnetic(3,3,nat)
  real(dp) :: sigma_paramagnetic(3,3,nat)
  real(dp) :: sigma_paramagnetic_us(3,3,nat)
  real(dp) :: sigma_paramagnetic_aug(3,3,nat)

  integer :: ik, ipol, jpol, i, ibnd, jbnd, isign, ispin
  real(dp) :: tmp(3,3), q(3), k_plus_q(3), braket, cc
  complex(dp), external :: zdotc

  !-----------------------------------------------------------------------
  ! allocate memory
  !-----------------------------------------------------------------------
  allocate ( p_evc(npwx,nbnd,3), vel_evc(npwx,nbnd,3) )
  allocate ( aux(npwx,nbnd), G_vel_evc(npwx,nbnd,3) )
  if (okvan) allocate ( svel_evc(npwx,nbnd,3), u_svel_evc(npwx,nbnd,3) )

  ! zero the f-sum rule
  f_sum(:,:) = 0.d0
  
  ! zero the Q tensors
  q_pGv(:,:,:) = 0.d0
  q_vGv(:,:,:) = 0.d0

  ! zero the current and the field
  j_bare(:,:,:,:) = (0.d0,0.d0)
  b_ind(:,:,:) = (0.d0,0.d0)
  
  ! zero the chemical shift
  sigma_bare = 0.d0
  sigma_diamagnetic = 0.d0
  sigma_paramagnetic = 0.d0
  sigma_paramagnetic_us = 0.d0
  sigma_paramagnetic_aug = 0.d0 

  write(stdout, '(5X,''Computing the magnetic susceptibility'',$)')
  write(stdout, '(5X,''isolve='',I1,4X,''ethr='',E10.4)') isolve, conv_threshold
  !====================================================================
  ! loop over k-points
  !====================================================================
  do ik = 1, nks
#ifdef __PARA
    if (me_pool == root_pool) &
    write(*, '(5X,''k-point #'',I5,'' of '',I5,6X,''pool #'',I3)') &
      ik, nks, my_pool_id+1
#else
    write(stdout, '(5X,''k-point #'',I5,'' of '',I5)') ik, nks
#endif
    current_k = ik
    current_spin = isk(ik)
    
    ! initialize at k-point k 
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(:) = g2kin(:) * tpiba2
    call init_us_2(npw, igk, xk(1,ik), vkb)
    
    ! read wfcs from file and compute becp
    call get_buffer (evc, nwordwfc, iunwfc, ik)

    ! this is the case q = 0
    q(:) = 0.d0

    if (job /= 'f-sum') call compute_u_kq(ik, q)
    call init_gipaw_2_no_phase (npw, igk, xk (1, ik), paw_vkb)
    call calbec (npw, paw_vkb, evc, paw_becp)

    ! compute the terms that do not depend on 'q':
    ! 1. the diamagnetic contribution to the field: Eq.(58) of [1]
    diamagnetic_corr_tensor = 0.0d0
    call diamagnetic_correction (diamagnetic_corr_tensor)
    sigma_diamagnetic = sigma_diamagnetic + diamagnetic_corr_tensor

    ! 2. the paramagnetic US augmentation: Eq.(30) of [2]
    if (okvan) then
      paramagnetic_corr_tensor_aug = 0.d0
      call paramagnetic_correction_aug (paramagnetic_corr_tensor_aug)
      sigma_paramagnetic_aug = sigma_paramagnetic_aug + paramagnetic_corr_tensor_aug
    endif

    ! compute p_k|evc>, v_{k,k}|evc>, G_k v_{k,k}|evc> and s_{k,k}|evc>
    call apply_operators
    if (okvan) then 
        evq(:,:) = evc(:,:)
        call apply_occ_occ_us
    endif

    !------------------------------------------------------------------
    ! f-sum rule
    !------------------------------------------------------------------
    do ipol = 1, 3 
      do jpol = 1, 3
        do ibnd = 1, nbnd_occ(ik)
          ! this is the "p-G-v" term
          braket = 2.d0*real(zdotc(npw, p_evc(1,ibnd,ipol), 1, G_vel_evc(1,ibnd,jpol), 1), dp)
          f_sum(ipol,jpol) = f_sum(ipol,jpol) + wg(ibnd,ik) * braket

          ! this is the "occ-occ" term
          if (okvan) then
              braket = zdotc(npw, p_evc(1,ibnd,ipol), 1, u_svel_evc(1,ibnd,jpol), 1)
              f_sum(ipol,jpol) = f_sum(ipol,jpol) + wg(ibnd,ik) * braket
              !!do jbnd = 1, nbnd_occ (ik)
              !!   braket = - zdotc(npw, evc(1,ibnd), 1, p_evc(1,jbnd,ipol), 1) * &
              !!              zdotc(npw, evc(1,jbnd), 1, svel_evc(1,ibnd,jpol), 1)
              !!enddo
         endif
        enddo
      enddo
    enddo

    !------------------------------------------------------------------
    ! pGv and vGv contribution to chi_{bare}
    !------------------------------------------------------------------
    if (job /= 'f-sum') then
      do i = 1, 3
        call add_to_tensor(q_pGv(:,:,0), p_evc, G_vel_evc)
        call add_to_tensor(q_vGv(:,:,0), vel_evc, G_vel_evc)
      enddo
    endif
    
    !------------------------------------------------------------------
    ! loop over -q and +q
    !------------------------------------------------------------------
    do isign = -1, 1, 2
      if (job == 'f-sum') cycle
      
      ! loop over cartesian directions
      do i = 1, 3
        ! set the q vector
        q(:) = 0.d0
        q(i) = dble(isign) * q_gipaw
        
        ! compute the wfcs at k+q
        call compute_u_kq(ik, q)
        
        ! compute p_k|evc>, v_k|evc> and G_{k+q} v_{k+q,k}|evc>
        call apply_operators
      
        k_plus_q(1:3) = xk(1:3,ik) + q(1:3)
        call init_gipaw_2_no_phase(npw, igk, k_plus_q, paw_vkb)

        ! pGv and vGv contribution to chi_bare
        call add_to_tensor(q_pGv(:,:,isign), p_evc, G_vel_evc)
        call add_to_tensor(q_vGv(:,:,isign), vel_evc, G_vel_evc)
        
        ! now the j_bare term 
        call add_to_current(j_bare(:,:,:,current_spin), evc, G_vel_evc)
        if (okvan) then
          call apply_occ_occ_us
          call add_to_current(j_bare(:,:,:,current_spin), evc, u_svel_evc)
        endif

        ! paramagnetic terms
        call paramagnetic_correction(paramagnetic_corr_tensor, paramagnetic_corr_tensor_us, &
             G_vel_evc, u_svel_evc, i)
        call add_to_sigma_para(paramagnetic_corr_tensor, sigma_paramagnetic)
        if (okvan) call add_to_sigma_para(paramagnetic_corr_tensor_us, sigma_paramagnetic_us)

      enddo  ! i=x,y,z
    enddo  ! isign

  enddo  ! ik
  
#ifdef __PARA
  ! reduce over G-vectors
  call mp_sum( f_sum, intra_pool_comm )
  call mp_sum( q_pGv, intra_pool_comm )
  call mp_sum( q_vGv, intra_pool_comm )
#endif
  
#ifdef __PARA
  ! reduce over k-points
  call mp_sum( f_sum, inter_pool_comm )
  call mp_sum( q_pGv, inter_pool_comm )
  call mp_sum( q_vGv, inter_pool_comm )
  call mp_sum( j_bare, inter_pool_comm )
  call mp_sum( sigma_diamagnetic, inter_pool_comm )
  call mp_sum( sigma_paramagnetic, inter_pool_comm )
  call mp_sum( sigma_paramagnetic_us, inter_pool_comm )
  call mp_sum( sigma_paramagnetic_aug, inter_pool_comm )
#endif
  
  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,'(5X,''End of magnetic susceptibility calculation'')')
  write(stdout,*)

  ! free memory as soon as possible
  deallocate( p_evc, vel_evc, aux, G_vel_evc )
  if (okvan) deallocate( svel_evc, u_svel_evc )
  
  ! f-sum rule
  if (iverbosity > 0) then
    write(stdout, '(5X,''f-sum rule:'')')
    write(stdout, tens_fmt) f_sum
  endif
  call symmatrix (f_sum)
  write(stdout, '(5X,''f-sum rule (symmetrized):'')')
  write(stdout, tens_fmt) f_sum
  if (job == 'f-sum') return

  ! F_{ij} = (2 - \delta_{ij}) Q_{ij}
  do ipol = 1, 3
    do jpol = 1, 3
      f_pGv(ipol,jpol,:) = 2.0_dp*q_pGv(ipol,jpol,:)
      if (ipol == jpol) f_pGv(ipol,jpol,:) = q_pGv(ipol,jpol,:)

      f_vGv(ipol,jpol,:) = 2.0_dp*q_vGv(ipol,jpol,:)
      if (ipol == jpol) f_vGv(ipol,jpol,:) = q_vGv(ipol,jpol,:)
    enddo
  enddo
  
  ! compute chi_bare both pGv and vGv terms
  chi_bare_pGv(:,:) = f_pGv(:,:,1) - 2.0_dp*f_pGv(:,:,0) + f_pGv(:,:,-1)
  chi_bare_pGv(:,:) = -0.5_dp * chi_bare_pGv(:,:) * alpha ** 2 &
       / ( q_gipaw * tpiba)**2
  if (iverbosity > 0) then
    write(stdout, '(5X,''chi_bare pGv (HH) in paratec units:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) / alpha ** 2
  endif
  call symmatrix (chi_bare_pGv)
  if (iverbosity > 0) then
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) / alpha ** 2
  endif
  
  chi_bare_vGv(:,:) = f_vGv(:,:,1) - 2.0_dp*f_vGv(:,:,0) + f_vGv(:,:,-1)
  chi_bare_vGv(:,:) = -0.5_dp * chi_bare_vGv(:,:) * alpha ** 2 &
       / ( q_gipaw * tpiba)**2
  if (iverbosity > 0) then
    write(stdout, '(5X,''chi_bare vGv (VV) in paratec units:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) / alpha ** 2
  endif
  call symmatrix(chi_bare_vGv)
  if (iverbosity > 0) then
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) / alpha ** 2
  endif

  ! convert from atomic units to 10^{-6} cm^3 / mol
  tmp(:,:) = chi_bare_pGv(:,:) * 1e6_dp * a0_to_cm**3.0_dp * avogadro
  write(stdout, '(5X,''chi_bare pGv (HH) in 10^{-6} cm^3/mol:'')')
  write(stdout, tens_fmt) tmp(:,:)

  tmp(:,:) = chi_bare_vGv(:,:) * 1e6_dp * a0_to_cm**3.0_dp * avogadro
  write(stdout, '(5X,''chi_bare vGv (VV) in 10^{-6} cm^3/mol:'')')
  write(stdout, tens_fmt) tmp(:,:)

  !--------------------------------------------------------------------
  ! now get the current, induced field and chemical shifts
  !--------------------------------------------------------------------
  chi_bare_pGv(:,:) = chi_bare_pGv(:,:) / omega
  j_bare(:,:,:,:) = j_bare(:,:,:,:) * alpha &
       / ( 2.0_dp * q_gipaw * tpiba * omega )

  ! either you symmetrize the current ...
  do ispin = 1, nspin
#ifdef __PARA
    call psymmetrize_field(j_bare(:,:,:,ispin), 1)
#else
    call symmetrize_field(j_bare(:,:,:,ispin), 1)
#endif
  enddo

  ! compute induced field
  do ipol = 1, 3
    call biot_savart(ipol)
  enddo

  ! write fields to disk
  do i = 1, nspin
    if (trim(filcurr) /= '') &
      call write_tensor_field(filcurr, i, j_bare(1,1,1,i))
  enddo
  if (trim(filfield) /= '') &
    call write_tensor_field(filfield, 0, b_ind_r)

  ! ... or you symmetrize the induced field
  !call symmetrize_field(b_ind_r,0)
  !call field_to_reciprocal_space

  if (job == 'nmr') then
    ! compute bare chemical shift and print all results
    call compute_sigma_bare( chi_bare_pGv, sigma_bare )
    call print_chemical_shifts(sigma_bare, sigma_diamagnetic, sigma_paramagnetic, &
                               sigma_paramagnetic_us, sigma_paramagnetic_aug)
  endif
  
  
CONTAINS



  !====================================================================
  ! compute p_k|evc>, v_{k+q,k}|evc>, G_{k+q} v_{k+q,k}|evc> and s_{k+q,k}|evc>
  !====================================================================
  SUBROUTINE apply_operators
    IMPLICIT NONE
    integer ipol

    p_evc(:,:,:) = (0.d0,0.d0)
    vel_evc(:,:,:) = (0.d0,0.d0)
    svel_evc(:,:,:) = (0.d0,0.d0)

    do ipol = 1, 3
      call apply_p(evc, p_evc(1,1,ipol), ik, ipol, q)
      call apply_vel(evc, vel_evc(1,1,ipol), ik, ipol, q)
      if (okvan) call apply_vel_NL('S', evc, svel_evc(1,1,ipol), ik, ipol, q)
      ! necessary because aux is overwritten by subroutine greenfunction
      aux(:,:) = vel_evc(:,:,ipol)
      call greenfunction(ik, aux, G_vel_evc(1,1,ipol), q)
    enddo
  END SUBROUTINE apply_operators



  !====================================================================
  ! compute |evq><evq|s_{k+q,k}|evc>
  !====================================================================
  SUBROUTINE apply_occ_occ_us
    IMPLICIT NONE
    integer ipol
    complex(dp), allocatable :: ps(:,:)

    allocate( ps(nbnd,nbnd) )
    do ipol = 1, 3
      ps = (0.d0,0.d0)
      aux(:,:) = svel_evc(:,:,ipol)
      CALL ZGEMM('C', 'N', nbnd_occ(ik), nbnd_occ(ik), npw, &
                (1.d0,0.d0), evq(1,1), npwx, aux(1,1), npwx, (0.d0,0.d0), &
                ps(1,1), nbnd)
#ifdef __PARA
      call mp_sum(ps, intra_pool_comm)
#endif
      aux = (0.d0,0.d0)
      CALL ZGEMM('N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
                (1.d0,0.d0), evq(1,1), npwx, ps(1,1), nbnd, (0.d0,0.d0), &
                aux(1,1), npwx)
      u_svel_evc(:,:,ipol) = -aux(:,:)
    enddo
    deallocate(ps)
  END SUBROUTINE apply_occ_occ_us



  !====================================================================
  ! add contribution the Q tensors
  ! Q_{\alpha,\beta} += <(e_i \times ul)_\alpha | (e_i \times ur)_\beta>
  !====================================================================
  SUBROUTINE add_to_tensor(qt, ul, ur)
    IMPLICIT NONE
    real(dp), intent(inout) :: qt(3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd,3), ur(npwx,nbnd,3)
    real(dp) :: braket
    integer :: ibnd, ia, ib, comp_ia, comp_ib, ind(3,3), mult(3,3)

    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)

    do ia = 1, 3    ! ia = alpha
      comp_ia = ind(ia,i)
      if (mult(ia,i) == 0) cycle

      do ib = 1, 3    ! ib = beta
        comp_ib = ind(ib,i)
        if (mult(ib,i) == 0) cycle

        do ibnd = 1, nbnd_occ(ik)
          braket = real(zdotc(npw, ul(1,ibnd,comp_ia), 1, &
                                   ur(1,ibnd,comp_ib), 1), dp)
          qt(ia,ib) = qt(ia,ib) + wg(ibnd,ik) * &
                      braket * mult(ia,i) * mult(ib,i)
        enddo  ! ibnd

      enddo  ! ib
    enddo  ! ia
  END SUBROUTINE add_to_tensor



  !====================================================================
  ! add contribution the the current
  ! j(r)_{\alpha,\beta} += <ul|J(r)|(B\times e_i \cdot ur)>
  !====================================================================
  SUBROUTINE add_to_current(j, ul, ur)
    IMPLICIT NONE
    real(dp), intent(inout) :: j(nrxxs,3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd), ur(npwx,nbnd,3)
    real(dp) :: braket, fact
    integer :: ibdir, icomp, ind(3,3), mult(3,3)

    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)

    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir, i)
      fact = real(mult(ibdir,i)*isign)
      call j_para(fact, ul(1,1), ur(1,1,icomp), ik, q, j(1,1,ibdir))
    enddo
  END SUBROUTINE add_to_current
  


  !====================================================================
  ! Add contribution to current: Eq.(46) of [1]
  !====================================================================
  SUBROUTINE add_to_sigma_para(paramagnetic_correction, sigma_paramagnetic)
    IMPLICIT NONE
    real(dp), intent(in) :: paramagnetic_correction(3,3,nat)
    real(dp), intent(inout) :: sigma_paramagnetic(3,3,nat)
    real(dp) :: fact
    integer :: ibdir, icomp, ipol, ind(3,3), mult(3,3)
    
    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)
    
    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir,i)
      fact = real(mult(ibdir,i)*isign)
      
      do ipol = 1, 3
         sigma_paramagnetic ( ipol, icomp, : ) &
              = sigma_paramagnetic ( ipol, icomp, : ) &
              + fact * paramagnetic_correction ( ipol, ibdir, : ) &
              / ( 2 * q_gipaw * tpiba )
      end do
    enddo
  END SUBROUTINE add_to_sigma_para


END SUBROUTINE suscept_crystal
