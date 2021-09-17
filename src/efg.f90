!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE efg
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the electric field gradient (EFG)
  !  
  USE kinds,                  ONLY : dp 
  USE io_global,              ONLY : stdout
  USE constants,              ONLY : pi, tpi, fpi, angstrom_au, rytoev, electronvolt_si
  USE fft_base,               ONLY : dffts
  USE ions_base,              ONLY : nat, atm, ityp, zv
  USE symme,                  ONLY : symtensor
  USE lsda_mod,               ONLY : nspin
  USE mp,                     ONLY : mp_sum
  USE gipaw_module,           ONLY : q_efg, iverbosity
  USE gipaw_results,          ONLY : res_efg

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE
  complex(dp), allocatable:: tmp(:,:,:)
  real(dp), allocatable:: aux(:), rho_s(:,:)
  real(dp), allocatable :: efg_bare(:,:,:), efg_ion(:,:,:)
  real(dp), allocatable :: zion(:), efg_gipaw(:,:,:), efg_tot(:,:,:)
  integer :: alpha, beta, na
  real(dp):: v(3), axis(3,3)
  real(dp) :: eta, Cq

  call start_clock('efg')

  allocate( efg_bare(3,3,nat), efg_ion(3,3,nat), zion(nat) )
  allocate( efg_gipaw(3,3,nat), efg_tot(3,3,nat) )

  ! calculate the bare contribution
  allocate( rho_s(dffts%nnr,2), aux(dffts%nnr) )
  call get_smooth_density(rho_s)
  aux(:) = rho_s(:,1)
  if (nspin == 2) aux(:) = aux(:) + rho_s(:,2)
  call efg_bare_el(aux, efg_bare)
  deallocate( rho_s, aux )

  ! calculate ionic contribution
  allocate ( tmp(nat,3,3) )
  call ewald_dipole (tmp, zv)
  do na = 1, nat
      efg_ion(:,:,na) = real(tmp(na,:,:), kind=dp)
  enddo
  deallocate( tmp )

  ! calculate GIPAW correction
  call get_rho_up_down
  call efg_correction(efg_gipaw)
  
  ! print results
  write(stdout,*)
  write(stdout,'(5X,''ELECTRIC FIELD GRADIENTS TENSORS IN Hartree/bohrradius^2:'')')
  write(stdout,*)

  if (iverbosity > 0) then  
      write(stdout,'(5X,''----- bare term -----'')')
      do na = 1, nat
        do beta = 1, 3
          write(stdout,1000) atm(ityp(na)), na, (efg_bare(alpha,beta,na), alpha=1,3)
        enddo
        write(stdout,*)
     enddo

      write(stdout,'(5X,''----- ionic term -----'')')
      do na = 1, nat
        do beta = 1, 3
          write(stdout,1000) atm(ityp(na)), na, (efg_ion(alpha,beta,na), alpha=1,3)
        enddo
        write(stdout,*)
     enddo

      write(stdout,'(5X,''----- GIPAW term -----'')')
      do na = 1, nat
        do beta = 1, 3
          write(stdout,1000) atm(ityp(na)), na, (efg_gipaw(alpha,beta,na), alpha=1,3)
        enddo
        write(stdout,*)
     enddo
  endif

  do na = 1, nat
     efg_tot(:,:,na) = efg_bare(:,:,na) + efg_ion(:,:,na) + efg_gipaw(:,:,na)
  enddo

  if (iverbosity > 0) then  
      write(stdout,'(5X,''----- total EFG -----'')')
      do na = 1, nat
        do beta = 1, 3
          write(stdout,1000) atm(ityp(na)), na, (efg_tot(alpha,beta,na), alpha=1,3)
        enddo
        write(stdout,*)
     enddo
  endif

  ! symmetrise efg_tensor
  call symtensor(nat, efg_tot)
  write(stdout,'(5X,''----- total EFG (symmetrized) -----'')')
  do na = 1, nat
    do beta = 1, 3
      write(stdout,1000) atm(ityp(na)), na, (efg_tot(alpha,beta,na), alpha=1,3)
    enddo
    write(stdout,*)
  enddo
1000 FORMAT(5X,A,I3,2X,3(F14.6,2X))

  ! output in magres format
  call output_magres_begin('efg')
  call output_magres_efg(efg_tot)
  call output_magres_end
  res_efg = efg_tot
  
  ! calculate principal axis and spectroscopic parameters
  write(stdout,*)
  write(stdout,'(5X,''NQR/NMR SPECTROSCOPIC PARAMETERS:'')')

  do na = 1, nat
      call principal_axis(efg_tot(:,:,na), v, axis)
      write(stdout,1001) atm(ityp(na)), na, 'Vxx=', v(1), 'axis=(', axis(1:3,1), ')'
      write(stdout,1001) atm(ityp(na)), na, 'Vyy=', v(2), 'axis=(', axis(1:3,2), ')'
      write(stdout,1001) atm(ityp(na)), na, 'Vzz=', v(3), 'axis=(', axis(1:3,3), ')'

      eta = 0.d0
      if (abs(v(3)) > 1d-5) eta = (v(1)-v(2)) / v(3)

      Cq = v(3) * q_efg(ityp(na)) * rytoev * 2.d0 * angstrom_au ** 2 &
           * electronvolt_si * 1.d18 / 6.62620d0
     
      write(stdout,1002) atm(ityp(na)), na, q_efg(ityp(na)), Cq, eta
      write(stdout,*)
  enddo

  call start_clock('efg')

1001 FORMAT(5X,A,I3,4X,A,F10.4,4X,A,3F10.6,A)
1002 FORMAT(5X,A,I3,2X,'Q=',F15.4,' 1e-30 m^2',2X,'Cq=',F15.4,' MHz',2X,'eta=',F8.5)

END SUBROUTINE efg



!-----------------------------------------------------------------------
SUBROUTINE get_smooth_density(rho)
  !-----------------------------------------------------------------------
  !
  ! ... Get the charge density on the smooth grid
  !  
  USE kinds,                  ONLY : dp 
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : inter_pool_comm
  USE lsda_mod,               ONLY : current_spin, isk, nspin
  USE wvfct,                  ONLY : nbnd, wg, g2kin, current_k
  USE gvecw,                  ONLY : gcutw
  USE klist,                  ONLY : nks, xk, igk_k, ngk
  USE gvect,                  ONLY : ngm, g
  USE wavefunctions,          ONLY : evc
  USE cell_base,              ONLY : tpiba2, omega
  USE io_files,               ONLY : nwordwfc, iunwfc
  USE buffers,                ONLY : get_buffer
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : invfft
  USE gipaw_module,           ONLY : job
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: rho(dffts%nnr,nspin)
  !-- local variables ----------------------------------------------------
  complex(dp) :: psic(dffts%nnr)
  integer :: ibnd, ik, is
  logical :: save_tg
  integer :: npw

  rho = (0.d0,0.d0)

  !!!! disable task groups in this routine, it's not time consuming
  !!!save_tg = dffts%have_task_groups
  !!!dffts%have_task_groups = .false.

  ! loop over k-points
  do ik = 1, nks
     current_k = ik
     current_spin = isk(ik)
     npw = ngk(ik)
 
     ! initialize at k-point k and read wfcs from file
     call gk_sort(xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
     call get_buffer (evc, nwordwfc, iunwfc, ik)

     ! loop over bands
     do ibnd = 1, nbnd
       psic(:) = (0.d0,0.d0)
       psic(dffts%nl(igk_k(1:npw,ik))) = evc(1:npw,ibnd)
       call invfft ('Wave', psic, dffts)
       rho(:,current_spin) = rho(:,current_spin) + wg(ibnd,ik) * &
                             (dble(psic(:))**2 + aimag(psic(:))**2) / omega
     enddo         
  enddo
#ifdef __MPI
  ! reduce over k-points
  call mp_sum( rho, inter_pool_comm )
#endif
  !!!dffts%have_task_groups = save_tg

  if (job == 'hyperfine' .or. job == 'mossbauer') then
    do is = 1, nspin
#ifdef __MPI
      call psymmetrize_rho_s(rho(1,is))
#else
      call symmetrize_rho_s(rho(1,is))
#endif
    enddo
  endif

END SUBROUTINE get_smooth_density



!-----------------------------------------------------------------------
SUBROUTINE efg_bare_el(rho, efg_bare)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the bare contribution to the EFG
  !  
  USE kinds,                  ONLY : dp 
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : intra_pool_comm
  USE constants,              ONLY : tpi, fpi
  USE gvecs,                  ONLY : ngms
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : fwfft
  USE gvect,                  ONLY : g, gg, gstart
  USE ions_base,              ONLY : nat, tau
  USE gipaw_module,           ONLY : job

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(in) :: rho(dffts%nnr)
  real(dp), intent(out) :: efg_bare(3,3,nat)  
  !-- local variables ----------------------------------------------------
  complex(dp), allocatable :: efg_g(:,:,:), rhoaux(:)
  integer :: alpha, beta, ig, na
  real(dp) :: arg, fac, e2, trace
  complex(dp) :: phase

  e2 = 1.0_dp  ! hartree
  fac = fpi * e2
  
  allocate( efg_g(ngms,3,3), rhoaux(dffts%nnr) )
  efg_g(:,:,:) = (0.d0,0.d0)
  rhoaux(:) = cmplx(rho(:), kind=dp)

  ! transform density to reciprocal space
  CALL fwfft ('Rho', rhoaux, dffts)
  
  ! electric field gradient in the G-space
  do ig = gstart, ngms
     trace = 1.d0/3.d0 * gg(ig)
     do alpha = 1, 3
        efg_g(ig,alpha,alpha) = -trace
        do beta = 1, 3
           efg_g(ig,alpha,beta) = ( efg_g(ig,alpha,beta) + &
                       g(alpha,ig)*g(beta,ig)) * fac * rhoaux(dffts%nl(ig)) / gg(ig)
        enddo
     enddo
  enddo
  
  ! fourier transform on the atomic position
  efg_bare(:,:,:) = 0.d0
  do alpha = 1, 3
     do beta = 1, 3
        do na = 1, nat
           do ig = gstart, ngms
              arg = sum(tau(1:3,na) * g(1:3,ig)) * tpi
              phase = cmplx(cos(arg),sin(arg), kind=dp)
              efg_bare(alpha,beta,na) = efg_bare(alpha,beta,na) + &
                   real(efg_g(ig,alpha,beta) * phase, kind=dp)
           enddo
        enddo
     enddo
  enddo
#ifdef __MPI
  call mp_sum( efg_bare, intra_pool_comm )
#endif

  ! opposite sign for hyperfine
  if (job == 'hyperfine') efg_bare(:,:,:) = -efg_bare(:,:,:)

  deallocate( efg_g, rhoaux )
  return
END SUBROUTINE efg_bare_el



!-----------------------------------------------------------------------
SUBROUTINE efg_correction(efg_corr_tens)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the GIPAW contribution to the EFG
  !  
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE kinds,                 ONLY : dp
  USE uspp,                  ONLY : ap
  USE parameters,            ONLY : ntypx
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : g, ngm
  USE klist,                 ONLY : nks, xk, igk_k, ngk
  USE cell_base,             ONLY : tpiba2
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                 ONLY : g2kin, current_k, wg
  USE gvecw,                 ONLY : gcutw
  USE lsda_mod,              ONLY : current_spin, isk
  USE wavefunctions,  ONLY : evc
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE becmod,                ONLY : calbec
  USE constants,             ONLY : pi, fpi
  USE buffers
  USE gipaw_module,          ONLY : job, nbnd_occ
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp,                    ONLY : mp_sum
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: efg_corr_tens(3,3,nat)
  !-- local variables ----------------------------------------------------
  integer :: j, nt, ibnd, il1, il2, ik, nbs1, nbs2, kkpsi
  integer :: lm, m1, m2, lm1, lm2, l1, l2, nrc
  integer :: ijkb0, ih, jh, na, ikb, jkb, r_first
  integer :: s_min, s_maj, s_weight
  complex(dp) :: bec_product 
  real(dp), allocatable :: at_efg(:,:,:), work(:)
  complex(dp), allocatable :: efg_corr(:,:)
  integer :: npw
  
  allocate( efg_corr(9,nat) )
  efg_corr = 0.0_dp

  allocate ( at_efg(paw_nkb,paw_nkb,ntypx) ) 
  at_efg = 0.0_dp
 
  ! Select majority and minority spin components
  call select_spin(s_min, s_maj)
  
  ! calculate radial integrals: <aephi|1/r^3|aephi> - <psphi|1/r^3|psphi>
  do nt = 1, ntyp
     kkpsi = paw_recon(nt)%aephi(1)%kkpsi
     allocate( work(kkpsi) )
     
     r_first = 1
     if ( abs ( rgrid(nt)%r(1) ) < 1d-8 ) r_first = 2
     
     do il1 = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(il1)%label%nrc

        do il2 = 1, paw_recon(nt)%paw_nbeta
           work = 0.0_dp

           do j = r_first, nrc
              work(j) = &
                   ( paw_recon(nt)%aephi(il1)%psi(j) &
                   * paw_recon(nt)%aephi(il2)%psi(j) &
                   - paw_recon(nt)%psphi(il1)%psi(j) &
                   * paw_recon(nt)%psphi(il2)%psi(j) ) &
                   / rgrid(nt)%r(j) ** 3
           enddo
           
           call simpson(nrc,work,rgrid(nt)%rab,at_efg(il1,il2,nt))
           
        enddo
     enddo
     
     deallocate ( work )
  enddo
  
  !  calculate the reconstruction part
  do ik = 1, nks
     current_k = ik
     current_spin = isk(ik)
     npw = ngk(ik)
 
     ! different sign for spins only in "hyperfine", not "efg"
     if ( current_spin == s_min .and. job == "hyperfine" ) then
        s_weight = -1
     else
        s_weight = +1
     endif
     
     call gk_sort ( xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin )
     call get_buffer ( evc, nwordwfc, iunwfc, ik)
     
     call init_gipaw_2 ( npw, igk_k(1,ik), xk(1,ik), paw_vkb )
     call calbec ( npw, paw_vkb, evc, paw_becp )
     
     do ibnd = 1, nbnd_occ(ik)
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
                       lm2 = m2 + l2**2 
                       
                       bec_product = paw_becp(jkb,ibnd) * conjg( paw_becp(ikb,ibnd) )

                       do lm = 5, 9
                          efg_corr(lm,na) = efg_corr(lm,na) &
                               + s_weight * bec_product &
                               * at_efg(nbs1,nbs2,nt) &
                               * ap(lm,lm1,lm2) * wg(ibnd,ik)
                       enddo
                    enddo
                 enddo
                 
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              endif
              
           enddo
        enddo
     enddo
     
  enddo
  
  ! For hyperfine interaction the function used is:
  !   (3r_ir_j/r^2 - delta_i,j) / r^3, for EFG the other way around
  if (job == "hyperfine") efg_corr = -efg_corr
  
  !  transform in cartesian coordinates
  efg_corr_tens(1,1,:) =  sqrt(3.0_dp) * efg_corr(8,:) - efg_corr(5,:)
  efg_corr_tens(2,2,:) = -sqrt(3.0_dp) * efg_corr(8,:) - efg_corr(5,:)
  efg_corr_tens(3,3,:) = 2.0_dp * efg_corr(5,:)
  efg_corr_tens(1,2,:) = sqrt(3.0_dp) * efg_corr(9,:)
  efg_corr_tens(2,1,:) = efg_corr_tens(1,2,:)
  efg_corr_tens(1,3,:) = -efg_corr(6,:) * sqrt(3.0_dp)
  efg_corr_tens(3,1,:) = efg_corr_tens(1,3,:)
  efg_corr_tens(2,3,:) = -efg_corr(7,:) * sqrt(3.0_dp)
  efg_corr_tens(3,2,:) = efg_corr_tens(2,3,:)
  
  efg_corr_tens = - sqrt(4.0_dp*pi/5.0_dp)*efg_corr_tens

#ifdef __MPI
  call mp_sum( efg_corr_tens, inter_pool_comm )
#endif
  
  deallocate ( efg_corr )
  deallocate ( at_efg )
  
END SUBROUTINE efg_correction

