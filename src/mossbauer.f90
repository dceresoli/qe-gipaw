!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE mossbauer
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the contact density at the nucleus.
  !  
  USE kinds,                  ONLY : dp 
  USE io_global,              ONLY : stdout
  USE parameters,             ONLY : ntypx
  USE fft_base,               ONLY : dffts
  USE ions_base,              ONLY : nat, atm, ityp, ntyp => nsp
  USE lsda_mod,               ONLY : nspin
 
  !-- local variables ----------------------------------------------------
  implicit none
  real(dp), allocatable :: rho_s(:,:), rho_tot(:)
  real(dp), allocatable :: moss_bare(:), moss_gipaw(:), moss_tot(:), moss_core(:)
  integer :: na

  call start_clock('mossbauer')

  !--------------------------------------------------------------------
  ! contact density
  !--------------------------------------------------------------------
  allocate( moss_bare(nat), moss_gipaw(nat), moss_tot(nat), moss_core(ntyp) )
  allocate( rho_s(dffts%nnr,nspin), rho_tot(dffts%nnr) )

  call get_smooth_density(rho_s)
  rho_tot(:) = rho_s(:,1)
  if (nspin == 2) rho_tot(:) = rho_tot(:) + rho_s(:,2)

  call moss_bare_el(rho_tot, moss_bare)
  deallocate( rho_s, rho_tot )

  call moss_gipaw_correction(moss_gipaw)
  moss_tot = moss_bare + moss_gipaw

  call moss_core_orbitals(moss_core)

  !--------------------------------------------------------------------
  ! Print results
  !--------------------------------------------------------------------
  write(stdout,*)
  write(stdout,'(5X,''VALENCE DENSITY AT NUCLEI (elec/bohrradius^3):'')')
  write(stdout,*)

  write(stdout,'(5X,8X,''  bare            core            GIPAW           total'')')
  do na = 1, nat
      write(stdout,1002) atm(ityp(na)), na, moss_bare(na), moss_core(ityp(na)), &
                         moss_gipaw(na), moss_tot(na) + moss_core(ityp(na))
  enddo
1002 FORMAT(5X,A,I3,2X,4(F14.6,2X))
  write(stdout,*)

  call stop_clock('mossbauer')
 
END SUBROUTINE mossbauer


!-----------------------------------------------------------------------
SUBROUTINE moss_core_orbitals(moss_core)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the core contribution to the Mossbauer density.
  ! ... For the time being, neglect relativistic and orbital contributions.
  !  
  USE kinds,                  ONLY : dp 
  USE constants,              ONLY : tpi, fpi
  USE ions_base,              ONLY : ntyp => nsp
  USE atom,                   ONLY : rgrid
  USE paw_gipaw,              ONLY : paw_recon

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: moss_core(ntyp)
  !-- local variables ----------------------------------------------------
  integer :: nt, core_orb, r_first, occupation
  real(dp) :: rho

  do nt = 1, ntyp
    moss_core(nt) = 0.0
    if ( paw_recon(nt)%gipaw_ncore_orbital == 0 ) cycle

    r_first = 1
    if ( abs ( rgrid(nt)%r(1) ) < 1d-8 ) r_first = 2

    do core_orb = 1, paw_recon(nt)%gipaw_ncore_orbital
      rho = paw_recon(nt)%gipaw_core_orbital(r_first,core_orb) ** 2 / rgrid(nt)%r(r_first)**2 / fpi
      occupation = 2 * ( 2 * paw_recon(nt)%gipaw_core_orbital_l(core_orb) + 1 )
      moss_core(nt) = moss_core(nt) + occupation * rho
    enddo

  enddo

END SUBROUTINE moss_core_orbitals


!-----------------------------------------------------------------------
SUBROUTINE moss_bare_el(rho_s, moss_bare)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the bare contribution to the Mossbauer density.
  ! ... For the time being, neglect relativistic and orbital contributions.
  !  
  USE kinds,                  ONLY : dp 
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : intra_pool_comm
  USE constants,              ONLY : tpi, fpi
  USE gvecs,                  ONLY : ngms
  USE gvect,                  ONLY : g, gstart
  USE parameters,             ONLY : ntypx
  USE ions_base,              ONLY : nat, tau
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : fwfft

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(in) :: rho_s(dffts%nnr)
  real(dp), intent(out) :: moss_bare(nat)
  !-- local variables ----------------------------------------------------
  complex(dp), allocatable :: rhoaux(:)
  integer :: ig, na
  real(dp) :: arg
  complex(dp) :: phase

  ! transform to reciprocal space
  allocate(rhoaux(dffts%nnr))
  rhoaux(:) = cmplx(rho_s(:), kind=dp)
  CALL fwfft('Rho', rhoaux, dffts)

  ! fourier transform on the atomic position
  moss_bare(:) = 0.d0
  do na = 1, nat
      do ig = gstart, ngms
          arg = sum(tau(1:3,na) * g(1:3,ig)) * tpi
          phase = cmplx(cos(arg),sin(arg), kind=dp)
          moss_bare(na) = moss_bare(na) + real(rhoaux(dffts%nl(ig)) * phase, kind=dp)
      enddo
  enddo
  call mp_sum(moss_bare, intra_pool_comm)

  deallocate(rhoaux)
  return
END SUBROUTINE moss_bare_el


!-----------------------------------------------------------------------
SUBROUTINE moss_gipaw_correction(moss_gipaw)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the GIPAW contribution to the Mossbauer contact density.
  ! ... For the time being, neglect relativistic and orbital contributions.
  !  
  USE kinds,                 ONLY : dp
  USE parameters,            ONLY : ntypx
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : g, ngm
  USE klist,                 ONLY : nks, xk, igk_k, ngk
  USE cell_base,             ONLY : tpiba2
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                 ONLY : g2kin
  USE gvecw,                 ONLY : gcutw
  USE wavefunctions,  ONLY : evc
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE becmod,                ONLY : calbec
  USE constants,             ONLY : pi, fpi
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp,                    ONLY : mp_sum
  USE buffers,               ONLY : get_buffer
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE lsda_mod,              ONLY : current_spin, isk
  USE wvfct,                 ONLY : current_k, wg
  USE gipaw_module,          ONLY : nbnd_occ, alpha
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: moss_gipaw(nat)
  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: at_moss(:,:,:)
  real(dp), allocatable :: work(:)
  integer :: j, nt, ibnd, il1, il2, ik, nbs1, nbs2, kkpsi
  integer :: m1, m2, lm1, lm2, l1, l2, nrc
  integer :: ijkb0, ih, jh, na, ikb, jkb
  integer :: r_first
  complex(dp) :: bec_product
  integer :: npw
  
  allocate( at_moss(paw_nkb,paw_nkb,ntyp) )
  at_moss = 0.0_dp
  
  ! calculate radial integrals: <aephi|aephi> - <psphi|psphi>
  do nt = 1, ntyp
     kkpsi = paw_recon(nt)%aephi(1)%kkpsi
     allocate( work(kkpsi) )

     r_first = 1
     if ( abs ( rgrid(nt)%r(1) ) < 1d-8 ) r_first = 2
     
     do il1 = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(il1)%label%nrc
        l1 = paw_recon(nt)%psphi(il1)%label%l
        if (l1 /= 0) cycle
        
        do il2 = 1, paw_recon(nt)%paw_nbeta
           l2 = paw_recon(nt)%psphi(il2)%label%l
           if (l2 /= 0) cycle
           
           work = 0.0_dp
           do j = r_first, nrc
                 work(j) = &
                      ( paw_recon(nt)%aephi(il1)%psi(j) &
                      * paw_recon(nt)%aephi(il2)%psi(j) &
                      - paw_recon(nt)%psphi(il1)%psi(j) &
                      * paw_recon(nt)%psphi(il2)%psi(j) ) &
                      / rgrid(nt)%r(j) ** 2 / fpi
           enddo
           
           ! density at the origin
           at_moss(il1,il2,nt) = work(r_first)
           
        enddo
     enddo

     deallocate ( work )
  enddo
  
  !  calculate the reconstruction part
  moss_gipaw = 0.d0
  
  do ik = 1, nks
     current_k = ik
     current_spin = isk(ik)
     npw = ngk(ik)
 
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
                    if (l1 /= 0) cycle
                    
                    do jh = 1, paw_recon(nt)%paw_nh
                       jkb = ijkb0 + jh
                       nbs2 = paw_recon(nt)%paw_indv(jh)
                       l2 = paw_recon(nt)%paw_nhtol(jh)
                       m2 = paw_recon(nt)%paw_nhtom(jh)
                       lm2 = m2 + l2**2 
                       if (l2 /= 0) cycle
                       
                       bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))
                       
                       moss_gipaw(na) = moss_gipaw(na) + at_moss(nbs1,nbs2,nt) * bec_product * wg(ibnd,ik)
                    enddo
                 enddo
                 
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              endif
              
           enddo
        enddo
     enddo
     
  enddo

  call mp_sum( moss_gipaw, inter_pool_comm )
  
END SUBROUTINE moss_gipaw_correction


