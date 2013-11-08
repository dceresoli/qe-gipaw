!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! NMR-related routines
!
! References:
!   [1] Phys. Rev. B 63, 245101 (2001)   (norm-conserving GIPAW)
!   [2] Phys. Rev. B 76, 024401 (2007)   (ultrasoft)
!   [3] Phys. Rev. B 76, 165122 (2007)   (metallic systems)


!====================================================================
! Paramagnetic contribution to the induced magnetic field
! norm-conserving contrib.: Eq.(44) of [1], Eq.(41) of [2] first term
! ultrasoft       contrib.: Eq.(41) of [2] second term
!====================================================================
SUBROUTINE paramagnetic_correction (paramagnetic_tensor, paramagnetic_tensor_us, &
                                    g_vel_evc, u_svel_evc, ipol)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, g2kin, current_k
  USE klist,                  ONLY : wk
  USE becmod,                 ONLY : calbec  
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp, paw_nkb, paw_recon
  USE gipaw_module,           ONLY : lx, ly, lz, paw_becp2, paw_becp3, alpha, &
                                     radial_integral_paramagnetic
  USE uspp,                   ONLY : okvan

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(inout) :: paramagnetic_tensor(3,3,nat)
  real(dp), intent(inout) :: paramagnetic_tensor_us(3,3,nat)
  complex(dp), intent(in) :: g_vel_evc(npwx,nbnd,3)
  complex(dp), intent(in) :: u_svel_evc(npwx,nbnd,3)
  integer, intent(in) :: ipol ! cartesian index of u_i

  !-- local variables ----------------------------------------------------
  complex(dp) :: para_corr(3,nat), para_corr_us(3,nat)
  complex(dp) :: bec_product, cc
  integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
  integer :: nt, ibnd, na, ijkb0, jpol

  do jpol = 1, 3 
     if ( jpol == ipol ) cycle
     call calbec (npw, paw_vkb, g_vel_evc(:,:,jpol), paw_becp2)
     if (okvan) call calbec (npw, paw_vkb, u_svel_evc(:,:,jpol), paw_becp3)

     para_corr = (0.d0, 0.d0)
     para_corr_us = (0.d0, 0.d0)

     do ibnd = 1, nbnd
        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat
              if (ityp (na) .eq.nt) then
                 do ih = 1, paw_recon(nt)%paw_nh
                    ikb = ijkb0 + ih
                    nbs1 = paw_recon(nt)%paw_indv(ih)
                    l1 = paw_recon(nt)%paw_nhtol(ih)
                    m1 = paw_recon(nt)%paw_nhtom(ih)
                    lm1 = m1 + l1**2
                      
                    do jh = 1, paw_recon(nt)%paw_nh
                       jkb = ijkb0 + jh
                       nbs2 = paw_recon(nt)%paw_indv(jh)
                       l2 = paw_recon(nt)%paw_nhtol(jh)
                       m2 = paw_recon(nt)%paw_nhtom(jh)
                       lm2 = m2 + l2**2
                         
                       if ( l1 /= l2 ) cycle
                       bec_product = conjg(paw_becp(ikb,ibnd)) * paw_becp2(jkb,ibnd)

                       cc = bec_product * radial_integral_paramagnetic(nbs1,nbs2,nt) &
                            !WAS: * wg(ibnd,current_k) * alpha**2
                            * wk(current_k) * alpha**2
                       para_corr(1,na) = para_corr(1,na) + cc * lx(lm1,lm2)
                       para_corr(2,na) = para_corr(2,na) + cc * ly(lm1,lm2)
                       para_corr(3,na) = para_corr(3,na) + cc * lz(lm1,lm2)

                       if (okvan) then
                         bec_product = conjg(paw_becp(ikb,ibnd)) * paw_becp3(jkb,ibnd)
                         cc = bec_product * radial_integral_paramagnetic(nbs1,nbs2,nt) &
                              !WAS:* wg (ibnd,current_k) * alpha ** 2
                              * wk(current_k) * alpha ** 2
                         para_corr_us(1,na) = para_corr_us(1,na) + cc * lx(lm1,lm2)
                         para_corr_us(2,na) = para_corr_us(2,na) + cc * ly(lm1,lm2)
                         para_corr_us(3,na) = para_corr_us(3,na) + cc * lz(lm1,lm2)
                       endif

                    enddo
                 enddo
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              endif
           enddo
        enddo
     enddo
     paramagnetic_tensor(1:3,jpol,1:nat) = real(para_corr(1:3,1:nat),dp)
     if (okvan) paramagnetic_tensor_us(1:3,jpol,1:nat) = real(para_corr_us(1:3,1:nat),dp)
  end do
END SUBROUTINE paramagnetic_correction


  
!====================================================================
! Diamagnetic contribution to the induced magnetic field
! Eq.(58) of [1], same as Eq.(43) of [2]
!====================================================================
SUBROUTINE diamagnetic_correction (diamagnetic_tensor)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k
  USE becmod,                 ONLY : calbec  
  USE constants,              ONLY : pi
  USE parameters,             ONLY : lmaxx
  USE uspp,                   ONLY : ap
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp, paw_nkb, paw_recon
  USE gipaw_module,           ONLY : lx, ly, lz, paw_becp2, alpha, &
                                     radial_integral_diamagnetic
  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(inout) :: diamagnetic_tensor(3,3,nat)

  !-- local variables ----------------------------------------------------
  integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
  integer :: nt, ibnd, na, lm, ijkb0
  complex(dp) :: dia_corr(lmaxx**2,nat)
  complex(dp) :: bec_product
  
  dia_corr = 0.0_dp
  do ibnd = 1, nbnd
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if ( ityp(na) == nt ) then
              do ih = 1, paw_recon(nt)%paw_nh
                 ikb = ijkb0 + ih
                 nbs1 = paw_recon(nt)%paw_indv(ih)
                 l1 = paw_recon(nt)%paw_nhtol(ih)
                 m1 = paw_recon(nt)%paw_nhtom(ih)
                 lm1 = m1 + l1**2
                 do jh = 1, paw_recon(nt)%paw_nh
                    jkb = ijkb0 + jh
                    nbs2 = paw_recon(nt)%paw_indv(jh)
                    l2 = paw_recon(nt)%paw_nhtol(jh)
                    m2 = paw_recon(nt)%paw_nhtom(jh)
                    lm2  = m2 + l2**2 
                      
                    bec_product = paw_becp(jkb,ibnd) * conjg( paw_becp(ikb,ibnd) )
                     
                    !<apsi> s/non-trace-zero component
                    ! 2/3 to separate the non-trace vanishing component
                    ! 1/(2c^2) from the equation (59) in PM-PRB
                    if ( l1 == l2 .and. m1 == m2 ) then
                       diamagnetic_tensor(1,1,na) &
                            = diamagnetic_tensor(1,1,na) &
                            + 2.0_dp / 3.0_dp * bec_product &
                            * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                            * wg(ibnd,current_k) * alpha ** 2 / 2.0_dp
                       diamagnetic_tensor(2,2,na) &
                            = diamagnetic_tensor(2,2,na) &
                            + 2.0_dp / 3.0_dp * bec_product &
                            * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                            * wg(ibnd,current_k) * alpha ** 2 / 2.0_dp
                       diamagnetic_tensor(3,3,na) &
                            = diamagnetic_tensor(3,3,na) &
                            + 2.0_dp / 3.0_dp * bec_product &
                            * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                            * wg(ibnd,current_k) * alpha ** 2 / 2.0_dp
                    endif
                      
                    ! 2/3 to separate the non-trace vanishing component
                    do lm = 5, 9
                       dia_corr(lm,na) =  dia_corr(lm,na) &
                            + bec_product / 3.0_dp &
                            * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                            * ap(lm,lm1,lm2) * wg(ibnd,current_k) * alpha ** 2 &
                            / 2.0_dp
                    enddo
                 enddo
              enddo
              ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
           endif
        enddo
     enddo
  enddo
    
  !  transform in cartesian coordinates
  dia_corr(5:9,:nat) = - sqrt(4.0_dp * pi/5.0_dp) * dia_corr(5:9,:nat)
  diamagnetic_tensor(1,1,:) = diamagnetic_tensor(1,1,:) &
       + sqrt(3.0_dp) * dia_corr(8,:) - dia_corr(5,:)
  diamagnetic_tensor(2,2,:) = diamagnetic_tensor(2,2,:) &
       - sqrt(3.0_dp) * dia_corr(8,:) - dia_corr(5,:)
  diamagnetic_tensor(3,3,:) = diamagnetic_tensor(3,3,:) &
       + dia_corr(5,:) * 2.0_dp
  diamagnetic_tensor(1,2,:) = diamagnetic_tensor(1,2,:) &
       +  dia_corr(9,:) * sqrt(3.0_dp)
  diamagnetic_tensor(2,1,:) = diamagnetic_tensor(1,2,:)
  diamagnetic_tensor(1,3,:) = diamagnetic_tensor(1,3,:) &
       - dia_corr(6,:) * sqrt(3.0_dp)
  diamagnetic_tensor(3,1,:) = diamagnetic_tensor(1,3,:)
  diamagnetic_tensor(2,3,:) = diamagnetic_tensor(2,3,:) &
       - dia_corr(7,:) * sqrt(3.0_dp)
  diamagnetic_tensor(3,2,:) = diamagnetic_tensor(2,3,:)
  ! dia_corr(5,:) = 3z^2-1
  ! dia_corr(6,:) = -xz
  ! dia_corr(7,:) = -yz
  ! dia_corr(8,:) = x^2-y^2
  ! dia_corr(9,:) = xy
END SUBROUTINE diamagnetic_correction



!====================================================================
! Ultrasoft augmentation (L_R Q_R) contribution to the bare and
! paramagnetic current
! TODO: modify for metals
!====================================================================
SUBROUTINE paramagnetic_correction_aug (paug_corr_tensor, j_bare_s)
  USE kinds,                  ONLY : dp
  USE ions_base,              ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                  ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                     current_k, ecutwfc
  USE lsda_mod,               ONLY : current_spin
  USE wavefunctions_module,   ONLY : evc
  USE becmod,                 ONLY : calbec, allocate_bec_type, deallocate_bec_type
  USE constants,              ONLY : pi
  USE parameters,             ONLY : lmaxx
  USE fft_base,               ONLY : dffts
  USE lsda_mod,               ONLY : nspin
  USE uspp,                   ONLY : ap
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp, paw_nkb, paw_recon
  USE gipaw_module,           ONLY : lx, ly, lz, radial_integral_paramagnetic, &
                                     q_gipaw, alpha, nbnd_occ, iverbosity
  USE uspp,                   ONLY : qq, vkb, nkb 
  USE uspp_param,             ONLY : nh
  USE cell_base,              ONLY : tpiba, omega, tpiba2
  USE klist,                  ONLY : xk
  USE gvect,                  ONLY : g, ngm
  USE io_global,              ONLY : stdout, ionode
#ifdef __BANDS
  USE gipaw_module,           ONLY : ibnd_start, ibnd_end
  USE mp,                     ONLY : mp_sum
  USE mp_bands,               ONLY : inter_bgrp_comm
#endif
  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(inout) :: paug_corr_tensor(3,3,nat)
  real(dp), intent(inout) :: j_bare_s(dffts%nnr,3,3,nspin)

  !-- local variables ----------------------------------------------------
  complex(dp), allocatable::pcorr_jpaug(:,:) 
  integer ::ibnd,ijkb0,nt,na,ih,ikb,jh,jkb,ig,ipol,jpol,kpol,kfac 
  integer :: nbs1,nbs2,l1,l2,m1,m2,lm1,lm2, ik
  complex(dp) , allocatable :: ps(:,:) !1st part
  complex(dp) , allocatable ::dvkbj(:,:),dvkby(:,:),Lp(:,:,:)!2nd part
  real(DP), allocatable  :: gk (:,:), gg(:,:)
  complex(dp) , allocatable :: LQ(:,:,:) !3rd part
  complex(dp) , allocatable :: aux1(:,:),aux2(:,:),paw_becp_gLQ(:,:)
  complex(dp) , allocatable :: g_LQ_evc(:,:,:),becp2(:,:)
  complex(dp) :: cc,bec_product
  real(dp) :: epsi(3,3), xyz(3,3),emine_q(3),dvkb_dir(3),ffact
  DATA epsi/0.d0,-3.d0,2.d0,3.d0,0.d0,-1.d0,-2.d0,1.d0,0.d0/, &
       xyz/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/

  !calculating ps = q_{ji}<p_i|u>
  emine_q(1)=0.d0; emine_q(2)=0.d0;emine_q(3)=0.d0

  allocate( ps( nkb, nbnd ), becp2(nkb,nbnd) )
  ik = current_k

  call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
  vkb = (0.d0,0.d0)
  call init_us_2 (npw, igk, xk(:,ik), vkb)
#ifdef __BANDS
  call calbec_bands (npwx, npw, nkb, vkb, evc, becp2, nbnd, ibnd_start, ibnd_end)
#else
  call calbec (npw, vkb, evc, becp2, nbnd)
#endif

  ps(:,:) = (0.d0, 0.d0)

  ijkb0 = 0 
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
#ifdef __BANDS
           do ibnd = ibnd_start, ibnd_end
#else
           do ibnd = 1, nbnd
#endif
              do ih = 1, nh(nt) 
                 ikb = ijkb0 + ih
                 do jh = 1, nh(nt)
                    jkb = ijkb0 + jh  
                    ps(jkb, ibnd) = ps(jkb, ibnd) + qq(jh,ih,nt) * becp2(ikb,ibnd)
                 enddo ! jh
              enddo ! ih
           enddo ! nbnd
           ijkb0 = ijkb0 + nh(nt)
        endif ! ityp(na)==nt
     enddo ! nat
  enddo ! ntyp

  ! now we have ps (nkb x nbnd)
  !calculating L|p> = eps(abc)r_b p_c |p> dir_a = eps(abc)p_c r_b|p>
  !                 = eps(abc)p_c (d|p>/dk_b) (derivative wrt k_b)
  !                 = eps(abc)(k+G)_c(d|p>/dk_b) (check the p part again later)
  allocate(dvkbj(npwx,nkb), dvkby(npwx,nkb), Lp(npwx,nkb,3))
  allocate(gk(3,npwx),gg (3,npwx))
  dvkbj(:,:) = (0.d0,0.d0); dvkby(:,:) = (0.d0,0.d0);Lp(:,:,:) = (0.d0,0.d0)
  gk(:,:) = 0.d0; gg(:,:) = 0.d0
  do ig = 1,npw
     gk(1:3,ig)=(xk(1:3,ik)+g(1:3,igk(ig)))*tpiba
     g2kin(ig) = SUM (gk(1:3,ig)**2)
     if(g2kin (ig) < 1.0d-10) then
       gg (:, ig) = 0.d0
     else
       gg (1:3,ig) = gk (1:3,ig) / sqrt(g2kin(ig))
     endif
  enddo
  !write(*,*)"ok 2.2"
  call gen_us_dj(ik,dvkbj)
  do ipol = 1,3
     dvkb_dir(:)= xyz(:,ipol)
     call gen_us_dy(ik,dvkb_dir, dvkby)
     do jpol = 1,3
        kpol = int(abs(epsi(ipol,jpol)))
        if(kpol.eq.0)cycle
        kfac = int(epsi(ipol,jpol)/kpol)
        ijkb0 =0
        do nt = 1,ntyp
           do na = 1,nat
              if (nt==ityp(na)) then
                 do ikb = 1, nh(nt)
                    ijkb0 = ijkb0+1
                    do ig = 1,npw
                       Lp(ig,ijkb0,kpol)=Lp(ig,ijkb0,kpol)+ &
                                (dvkby(ig,ijkb0) + dvkbj(ig,ijkb0) &
                                 *gg(ipol,ig))*gk(jpol,ig)*kfac
                    enddo!npw
                 enddo !ikb
              endif !ityp(na)=nt
           enddo !na
        enddo !nt
     enddo !jpol
  enddo !ipol 
  ! now we have both ps and Lp (npwx x nkb,3)
  ! we can construct LQ = LQ|u> = L|p>q<p|u>
  
  allocate (LQ(npwx,nbnd,3))
  LQ(:,:,:) = (0.d0,0.d0)
  do kpol = 1,3
#ifdef __BANDS
        call zgemm ('N', 'N', npwx, ibnd_end-ibnd_start+1, nkb, &
            (1.d0,0.d0),Lp(:,:,kpol), npwx, ps(1,ibnd_start), nkb, (1.d0,0.d0), &
            LQ(1,ibnd_start,kpol), npwx )
#else
        call zgemm ('N', 'N', npwx, nbnd, nkb, &
            (1.d0,0.d0),Lp(:,:,kpol), npwx, ps, nkb, (1.d0,0.d0), &
            LQ(1,1,kpol), npwx )
#endif
  enddo
#if defined(__MPI) && defined(__BANDS) 
  call mp_sum(LQ, inter_bgrp_comm)
#endif
  ! now we have LQ (npw,nbnd)  
  !apply Greens function
  allocate(aux1(npwx,nbnd),aux2 (npwx,nbnd),pcorr_jpaug(3,nat))
  allocate(g_LQ_evc(npwx,nbnd,3),paw_becp_gLQ(paw_nkb,nbnd))
  do kpol = 1,3
     pcorr_jpaug(:,:) = (0.0d0,0.d0)
     aux1(:,:) = (0.d0,0.d0);aux2(:,:) = (0.d0,0.d0)
     aux1(:,:) = LQ(:,:,kpol)!coz it changes in gf
 call start_clock('para:gf')
     call greenfunction(ik, aux1, aux2, (/0.d0, 0.d0, 0.d0/))
 call stop_clock('para:gf')
     g_LQ_evc(:,:,kpol) = aux2(:,:)
#ifdef __BANDS
     call calbec_bands (npwx, npw, paw_nkb, paw_vkb, aux2, paw_becp_gLQ, nbnd,ibnd_start, ibnd_end)
     do ibnd = ibnd_start, ibnd_end
#else
     call calbec (npw, paw_vkb , aux2 ,paw_becp_gLQ)
     do ibnd = 1, nbnd
#endif
        ijkb0 = 0
        do nt = 1 , ntyp
              do na = 1, nat
                 if(ityp(na).eq.nt) then
                   do ih = 1, paw_recon(nt)%paw_nh
                      ikb = ijkb0 + ih
                      nbs1 = paw_recon(nt)%paw_indv(ih);l1 = paw_recon(nt)%paw_nhtol(ih)
                      m1 = paw_recon(nt)%paw_nhtom(ih); lm1 = m1 + l1**2
                      do jh = 1, paw_recon(nt)%paw_nh
                         jkb = ijkb0 + jh
                         nbs2 = paw_recon(nt)%paw_indv(jh);l2 = paw_recon(nt)%paw_nhtol(jh)
                         m2 = paw_recon(nt)%paw_nhtom(jh);lm2=m2+l2**2
                         if(l1 /= l2) cycle !not sure of this..
                         bec_product = conjg(paw_becp(ikb,ibnd))*paw_becp_gLQ(jkb,ibnd)
                         cc = bec_product*radial_integral_paramagnetic(nbs1,nbs2,nt) &
                              * wg(ibnd,ik)* alpha ** 2
                         pcorr_jpaug(1,na) = pcorr_jpaug(1,na) + cc * lx ( lm1, lm2 )
                         pcorr_jpaug(2,na) = pcorr_jpaug(2,na) + cc * ly ( lm1, lm2 )
                         pcorr_jpaug(3,na) = pcorr_jpaug(3,na) + cc * lz ( lm1, lm2 )
                      enddo !jh
                   enddo !ih
                   ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
                 endif !ityp(na)==nt
             enddo ! nat
        enddo !ntyp 
     enddo !bands
#if defined(__MPI) && defined(__BANDS)
     call mp_sum(pcorr_jpaug, inter_bgrp_comm)
#endif  
     paug_corr_tensor(:,kpol,:) = REAL (pcorr_jpaug(:,:), dp)
     if ( iverbosity > 20 ) then
        write(6,'("PARA_AUG",1I3,3(F16.7,2X))') &
             kpol,  REAL ( pcorr_jpaug(1:3,1) ) * 1e6
     endif
  enddo !kpol   
    

  ffact =  -( 2.0_dp * q_gipaw * tpiba )
  do ih = 1,3
        call j_para(ffact,evc,g_LQ_evc(:,:,ih),ik,emine_q,j_bare_s(:,:,ih,current_spin))
  enddo
  ! adding this to sigma is easy coz there is no cross product..
  deallocate(pcorr_jpaug)
END SUBROUTINE paramagnetic_correction_aug

  

!====================================================================
! Compute the bare contribution to the chemical shift by evaluating
! the induced magnetic field at the nuclei
!====================================================================
SUBROUTINE compute_sigma_bare(B_ind, chi_bare, sigma_bare, sigma_shape)
  USE kinds,                ONLY : dp
  USE gvect,                ONLY : ngm, gstart, nl, nlm, g
  USE ions_base,            ONLY : nat, tau, atm, ityp
  USE pwcom,                ONLY : pi, tpi
  USE gipaw_module,         ONLY : use_nmr_macroscopic_shape, &
                                   nmr_macroscopic_shape
#ifdef __BANDS
  USE mp_bands,             ONLY : intra_bgrp_comm
#else
  USE mp_pools,             ONLY : intra_pool_comm
#endif
  USE lsda_mod,             ONLY : nspin
  USE mp,                   ONLY : mp_sum
  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  complex(dp), intent(in) :: B_ind(ngm,3,3,nspin)
  real(dp), intent(in) :: chi_bare(3,3)
  real(dp), intent(out) :: sigma_bare(3,3,nat), sigma_shape(3,3)
  !-- local variables ---------------------------------------------------
  integer :: na, ig, ispin
  real(dp) :: arg
  complex(dp) :: tmp_sigma(3,3)
  
  do na = 1, nat
    tmp_sigma(:,:) = 0.d0
    
    do ispin = 1, nspin
      do ig = gstart, ngm
        arg = (g(1,ig)*tau(1,na) + g(2,ig)*tau(2,na) + g(3,ig)*tau(3,na)) * tpi
        tmp_sigma(:,:) = tmp_sigma(:,:) + B_ind(ig,:,:,ispin) * cmplx(cos(arg),sin(arg), dp)
      enddo
    enddo
   
    sigma_bare(:,:,na) = real(tmp_sigma(:,:), dp)
  enddo
#ifdef __MPI
#ifdef __BANDS
  call mp_sum( sigma_bare, intra_bgrp_comm )
#else
  call mp_sum( sigma_bare, intra_pool_comm )
#endif
#endif

  if (use_nmr_macroscopic_shape) then
    sigma_shape(:,:) = -4.d0*pi * nmr_macroscopic_shape(:,:) * chi_bare(:,:)
  else
    sigma_shape(:,:) = 0.d0
  end if

END SUBROUTINE compute_sigma_bare



!====================================================================
! Print the contributions to the chemical shift and the total sigma
!====================================================================
SUBROUTINE print_chemical_shifts(sigma_shape, sigma_bare, sigma_diamagnetic, sigma_paramagnetic, &
                                 sigma_paramagnetic_us, sigma_paramagnetic_aug, sigma_tot)
  USE kinds,                ONLY : dp
  USE ions_base,            ONLY : nat, tau, atm, ityp, ntyp => nsp
  USE io_global,            ONLY : stdout
  USE symme,                ONLY : symtensor
  USE uspp,                 ONLY : okvan
  USE pwcom,                ONLY : lgauss
  USE gipaw_module,         ONLY : tens_fmt, iverbosity, nmr_shift_core
  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(in) :: sigma_shape(3,3), sigma_bare(3,3,nat)
  real(dp), intent(inout) :: sigma_diamagnetic(3,3,nat)
  real(dp), intent(inout) :: sigma_paramagnetic(3,3,nat)
  real(dp), intent(inout) :: sigma_paramagnetic_us(3,3,nat)
  real(dp), intent(inout) :: sigma_paramagnetic_aug(3,3,nat)
  real(dp), intent(out) :: sigma_tot(3,3,nat)
  !-- local variables ---------------------------------------------------
  real(dp) :: tr_sigma, axis(3,3), v(3), aniso, eta
  integer :: na
  
  write(stdout,'(5X,''Contributions to the NMR chemical shifts: -------------------------------'')')
  write(stdout,*)

  call trace(3, sigma_shape, tr_sigma)
  tr_sigma = tr_sigma / 3.d0
  write(stdout,'(5X,''Macroscopic shape contribution in ppm:'',10X,F14.2)') tr_sigma*1.0d6
  if (iverbosity > 0) write(stdout,tens_fmt) sigma_shape(:,:) * 1.0d6
  write(stdout,*)

  write(stdout,'(5X,''Core contribution in ppm:'')')
  write(stdout,*)
  do na = 1, nat
    write(stdout,'(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  core sigma: '',F14.2)') &
          na, atm(ityp(na)), tau(:,na), nmr_shift_core(ityp(na))*1.0d6
  enddo
  write(stdout,*)

  write(stdout,'(5X,''Bare contribution in ppm:'')')
  write(stdout,*)
  do na = 1, nat
    call trace(3, sigma_bare(1,1,na), tr_sigma)
    tr_sigma = tr_sigma / 3.d0
    write(stdout,'(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  bare sigma: '',F14.2)') &
          na, atm(ityp(na)), tau(:,na), tr_sigma * 1.0d6
    if (iverbosity > 0) write(stdout,tens_fmt) sigma_bare(:,:,na) * 1.0d6
  enddo
  
  write(stdout,'(5X,''Diamagnetic contribution in ppm:'')')
  write(stdout,*)
  call symtensor (nat, sigma_diamagnetic)
  do na = 1, nat
    call trace(3, sigma_diamagnetic(1,1,na), tr_sigma)
    tr_sigma = tr_sigma / 3.d0
    write(stdout,'(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  dia sigma: '',F14.2)') &
        na, atm(ityp(na)), tau(:,na), tr_sigma*1.0d6
    if (iverbosity > 0) write(stdout,tens_fmt) sigma_diamagnetic(:,:,na) * 1.0d6
  enddo

  write(stdout,'(5X,''Paramagnetic contribution in ppm:'')')
  write(stdout,*)
  call symtensor (nat, sigma_paramagnetic)
  do na = 1, nat
    call trace(3, sigma_paramagnetic(1,1,na), tr_sigma)
    tr_sigma = tr_sigma / 3.d0
    write(stdout,'(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  para sigma: '',F14.2)') &
        na, atm(ityp(na)), tau(:,na), tr_sigma*1.0d6
    if (iverbosity > 0) write(stdout,tens_fmt) sigma_paramagnetic(:,:,na) * 1.0d6
  enddo

  if (okvan) then
    write(stdout,'(5X,''Paramagnetic US occ-occ contribution in ppm:'')')
    write(stdout,*)
    call symtensor (nat, sigma_paramagnetic_us)
    do na = 1, nat
      call trace(3, sigma_paramagnetic_us(1,1,na), tr_sigma)
      tr_sigma = tr_sigma / 3.d0
      write(stdout,'(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  para_oo sigma: '',F14.2)') &
          na, atm(ityp(na)), tau(:,na), tr_sigma*1.0d6
      if (iverbosity > 0) write(stdout,tens_fmt) sigma_paramagnetic_us(:,:,na) * 1.0d6
    enddo

    write(stdout,'(5X,''Paramagnetic US L_R Q_R contribution in ppm:'')')
    write(stdout,*)
    call symtensor (nat, sigma_paramagnetic_aug)
    do na = 1, nat
      call trace(3, sigma_paramagnetic_aug(1,1,na), tr_sigma)
      tr_sigma = tr_sigma / 3.d0
      write(stdout,'(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  para_lq sigma: '',F14.2)') &
          na, atm(ityp(na)), tau(:,na), tr_sigma*1.0d6
      if (iverbosity > 0) write(stdout,tens_fmt) sigma_paramagnetic_aug(:,:,na) * 1.0d6
    enddo
  endif

  write(stdout,*)
  write(stdout,'(5X,''Total NMR chemical shifts in ppm: ---------------------------------------'')')
  write(stdout,'(5X,''(adopting the Simpson convention for anisotropy and asymmetry)-----------'')')
  write(stdout,*)
  sigma_tot = sigma_bare + sigma_diamagnetic + sigma_paramagnetic
  if (okvan) sigma_tot = sigma_tot + sigma_paramagnetic_us + sigma_paramagnetic_aug
  do na = 1, nat
    sigma_tot(1,1,na) = sigma_tot(1,1,na) + nmr_shift_core(ityp(na))
    sigma_tot(2,2,na) = sigma_tot(2,2,na) + nmr_shift_core(ityp(na))
    sigma_tot(3,3,na) = sigma_tot(3,3,na) + nmr_shift_core(ityp(na))
    sigma_tot(:,:,na) = sigma_tot(:,:,na) + sigma_shape

    call trace(3, sigma_tot(1,1,na), tr_sigma)
    tr_sigma = tr_sigma / 3.d0
    write(stdout,'(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  Total sigma: '',F14.2)') &
        na, atm(ityp(na)), tau(:,na), tr_sigma*1.0d6
    if (iverbosity > 0) write(stdout,tens_fmt) sigma_tot(:,:,na) * 1.0d6

    call principal_axis_simpson(sigma_tot(:,:,na), v, axis)
    aniso = v(3) - tr_sigma
    if (abs(aniso) > 1d-6) then
      eta = (v(2) - v(1))/aniso
    else
      eta = 0.d0
    endif
    aniso = aniso * 1.5d0
    write(stdout,1000) atm(ityp(na)), na, 'anisotropy:', aniso*1.0d6, 'eta:', eta
    write(stdout,1001) atm(ityp(na)), na, 'sigma_11=', v(1)*1.0d6, 'axis=(', axis(1:3,1), ')'
    write(stdout,1001) atm(ityp(na)), na, 'sigma_22=', v(2)*1.0d6, 'axis=(', axis(1:3,2), ')'
    write(stdout,1001) atm(ityp(na)), na, 'sigma_33=', v(3)*1.0d6, 'axis=(', axis(1:3,3), ')'
    write(stdout,*)
  enddo

  if (lgauss) &
    write(stdout,'(5X,''*** ATTENTION: system is metallic, Knight shift not included ***'')')

1000 FORMAT(5X,A,I3,4X,A,F10.2,4X,A,F10.4)
1001 FORMAT(5X,A,I3,4X,A,F10.4,4X,A,3F10.6,A)

END SUBROUTINE print_chemical_shifts

