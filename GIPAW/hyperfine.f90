!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! for the time being, neglect ZORA for heavy nuclei
#undef ZORA

!-----------------------------------------------------------------------
SUBROUTINE hyperfine
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the nuclear hyperfine couplings
  !  
  USE kinds,        ONLY : dp 
  USE io_global,    ONLY : stdout
  USE parameters,   ONLY : ntypx
  USE constants,    ONLY : pi, tpi, fpi, angstrom_au, rytoev, electronvolt_si, c_si
  USE gsmooth,      ONLY : nrxxs
  USE scf,          ONLY : rho
  USE symme,        ONLY : symtensor
  USE lsda_mod,     ONLY : current_spin, nspin
  USE wvfct,        ONLY : current_k
  USE ions_base,    ONLY : nat, tau, atm, ityp
  use constants,    ONLY : bohr_radius_si
  USE mp_global,    ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE gipaw_module, ONLY : hfi_nuclear_g_factor, hfi_output_unit, &
                           hfi_isotope, job, iverbosity, radial_integral_splines
 
  !-- constants ----------------------------------------------------------
  IMPLICIT NONE
  real(dp), parameter :: mu0_by_fpi = 1e-7
  real(dp), parameter :: mu_n = 5.05078324e-27_dp
  real(dp), parameter :: bohr_radius = bohr_radius_si
  real(dp), parameter :: gamma_e = 28024953.64_dp
  real(dp), parameter :: lambda = C_SI / 1.0e+8_dp
  real(dp), parameter :: common_factor = mu0_by_fpi * mu_n / Bohr_radius ** 3

  !-- local variables ----------------------------------------------------
  complex(dp), allocatable :: spin_den(:), rho_s(:,:)
  real(dp), allocatable :: hfi_dip_bare(:,:,:), hfi_dip_gipaw(:,:,:)
  real(dp), allocatable :: hfi_dip_tot(:,:,:)
  real(dp), allocatable :: hfi_fc_bare(:), hfi_fc_bare_zora(:)
  real(dp), allocatable :: hfi_fc_gipaw(:), hfi_fc_gipaw_zora(:)
  real(dp), allocatable :: hfi_fc_core(:), hfi_fc_tot(:)
  integer :: s_min, s_maj
  integer :: na, alpha, beta
  real(dp) :: output_factor, fact
  real(dp):: v(3), axis(3,3)

  call start_clock('hyperfine')

  ! select majority and minority spin components
  call select_spin(s_min, s_maj)

  !--------------------------------------------------------------------
  ! Dipolar contribution to hyperfine
  !--------------------------------------------------------------------
  allocate( hfi_dip_bare(3,3,nat), hfi_dip_gipaw(3,3,nat), hfi_dip_tot(3,3,nat) )

  ! calculate the bare dipole contribution
  allocate( rho_s(nrxxs,2), spin_den(nrxxs) )
  call get_smooth_density(rho_s)  ! this subroutine is efg.g90
  spin_den(:) = rho_s(:,s_maj) - rho_s(:,s_min)
  call efg_bare_el(spin_den, hfi_dip_bare)
  deallocate( rho_s)

  ! calculate GIPAW dipole correction
  call efg_correction(hfi_dip_gipaw)


  !--------------------------------------------------------------------
  ! Fermi-contact contribution to hyperfine
  !--------------------------------------------------------------------
  allocate( hfi_fc_bare(nat), hfi_fc_bare_zora(nat) )
  allocate( hfi_fc_gipaw(nat), hfi_fc_gipaw_zora(nat) )
  allocate( hfi_fc_core(nat), hfi_fc_tot(nat) )

  ! calculate the bare Fermi-contact contribution
  spin_den(:) = rho%of_r(:,s_maj) - rho%of_r(:,s_min)
  call hfi_fc_bare_el(spin_den, hfi_fc_bare, hfi_fc_bare_zora)
  deallocate( spin_den )

  ! calculate the GIPAW Fermi-contact contribution
  call hfi_fc_gipaw_correction(hfi_fc_gipaw, hfi_fc_gipaw_zora)

  ! calculate the core-relaxation Fermi-contact contribution
  call hfi_fc_core_relax(hfi_fc_core)


  !--------------------------------------------------------------------
  ! Print results
  !--------------------------------------------------------------------
  ! select the output unit
  select case (trim(hfi_output_unit))
      case('MHz')
          output_factor = 1d-3 * common_factor * gamma_e
      case('mT')
          output_factor = 1d3  * common_factor
      case('G', 'Gauss')
          output_factor = 1d4  * common_factor
      case('10e-4cm^-1')
          output_factor = 1d-3 * common_factor * gamma_e / lambda
      case default
          call errore( "hyperfine", "unknown units for output", 1)
  end select

  write(stdout,*)
  write(stdout,'(5X,''NUCLEAR G-TENSORS FROM INPUT:'')')
  do na = 1, nat
      write(stdout,'(5X,A,I3,2X,F12.6)') atm(ityp(na)), na, hfi_nuclear_g_factor(ityp(na))
  enddo

  write(stdout,*)
  write(stdout,'(5X,''NUCLEAR HYPERFINE COUPLINGS IN '',A,'':'')') hfi_output_unit
  write(stdout,*)

  write(stdout,'(5X,''DIPOLAR COUPLINGS:'')')
  if (iverbosity > 0) then  
      write(stdout,'(5X,''----- bare term -----'')')
      do na = 1, nat
        do beta = 1, 3
          fact = hfi_nuclear_g_factor(ityp(na)) * output_factor
          write(stdout,1000) atm(ityp(na)), na, (fact * hfi_dip_bare(alpha,beta,na), alpha=1,3)
        enddo
        write(stdout,*)
      enddo

      write(stdout,'(5X,''----- GIPAW term -----'')')
      do na = 1, nat
        do beta = 1, 3
          fact = hfi_nuclear_g_factor(ityp(na)) * output_factor
          write(stdout,1000) atm(ityp(na)), na, (fact * hfi_dip_gipaw(alpha,beta,na), alpha=1,3)
        enddo
        write(stdout,*)
      enddo
  endif

  do na = 1, nat
     hfi_dip_tot(:,:,na) = hfi_dip_bare(:,:,na) + hfi_dip_gipaw(:,:,na)
  enddo

  if (iverbosity > 0) then  
      write(stdout,'(5X,''----- total dipolar -----'')')
      do na = 1, nat
        do beta = 1, 3
          fact = hfi_nuclear_g_factor(ityp(na)) * output_factor
          write(stdout,1000) atm(ityp(na)), na, (fact * hfi_dip_tot(alpha,beta,na), alpha=1,3)
        enddo
        write(stdout,*)
      enddo
  endif

  call symtensor(nat, hfi_dip_tot)
  write(stdout,'(5X,''----- total dipolar (symmetrized) -----'')')
  do na = 1, nat
    do beta = 1, 3
      fact = hfi_nuclear_g_factor(ityp(na)) * output_factor
      write(stdout,1000) atm(ityp(na)), na, (fact * hfi_dip_tot(alpha,beta,na), alpha=1,3)
    enddo
    write(stdout,*)
  enddo
1000 FORMAT(5X,A,I3,2X,3(F14.6,2X))

  write(stdout,*)
  write(stdout,'(5X,''PRINCIPAL AXIS OF THE DIPOLAR COUPLINGS:'')')
  do na = 1, nat
      fact = hfi_nuclear_g_factor(ityp(na)) * output_factor
      call principal_axis(hfi_dip_tot(:,:,na)*fact, v, axis)
      write(stdout,1001) atm(ityp(na)), na, 'Axx=', v(1), 'axis=(', axis(1:3,1), ')'
      write(stdout,1001) atm(ityp(na)), na, 'Ayy=', v(2), 'axis=(', axis(1:3,2), ')'
      write(stdout,1001) atm(ityp(na)), na, 'Azz=', v(3), 'axis=(', axis(1:3,3), ')'
      write(stdout,*)
  enddo
1001 FORMAT(5X,A,I3,4X,A,F10.4,4X,A,3F10.6,A)


  if (iverbosity > 1) then
    write(stdout,*)
    write(stdout,'(5X,''SPIN DENSITIES IN bohrradius^-3 WITHOUT ZORA:'')')
    write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
    write(stdout,'(5X,8X,''  bare            GIPAW           core-relax      total'')')
    do na = 1, nat
        hfi_fc_tot(na) = hfi_fc_bare(na) + hfi_fc_gipaw(na) + hfi_fc_core(na)
        write(stdout,1002) atm(ityp(na)), na, hfi_fc_bare(na), &
            hfi_fc_gipaw(na), hfi_fc_core(na), hfi_fc_tot(na)
    enddo
  endif

  write(stdout,*)
  write(stdout,'(5X,''ISOTROPIC (FERMI-CONTACT) COUPLINGS WITHOUT ZORA:'')')
  write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
  write(stdout,'(5X,8X,''  bare            GIPAW           core-relax      total'')')
  do na = 1, nat
      hfi_fc_tot(na) = hfi_fc_bare(na) + hfi_fc_gipaw(na) + hfi_fc_core(na)
      fact = 8*pi/3 * output_factor * hfi_nuclear_g_factor(ityp(na))
      write(stdout,1002) atm(ityp(na)), na, fact * hfi_fc_bare(na), &
          fact * hfi_fc_gipaw(na), fact * hfi_fc_core(na), fact * hfi_fc_tot(na)
  enddo
1002 FORMAT(5X,A,I3,2X,4(F14.6,2X))

#ifdef ZORA
  write(stdout,*)
  write(stdout,'(5X,''ISOTROPIC (FERMI-CONTACT) COUPLINGS INCLUDING ZORA:'')')
  write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
  write(stdout,'(5X,8X,''  bare            GIPAW           core-relax      total'')')
  do na = 1, nat
      hfi_fc_tot(na) = hfi_fc_bare_zora(na) + hfi_fc_gipaw_zora(na) + hfi_fc_core(na)
      fact = 8*pi/3 * output_factor * hfi_nuclear_g_factor(ityp(na))
      write(stdout,1002) atm(ityp(na)), na, fact * hfi_fc_bare_zora(na), &
          fact * hfi_fc_gipaw_zora(na), fact * hfi_fc_core(na), fact * hfi_fc_tot(na)
  enddo
#endif

  call stop_clock('hyperfine')
 
END SUBROUTINE hyperfine



!-----------------------------------------------------------------------
SUBROUTINE hfi_fc_bare_el(rho_s, hfi_bare, hfi_bare_zora)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the bare contribution to the hyperfine
  !  
  USE kinds,        ONLY : dp 
  USE mp,           ONLY : mp_sum
  USE mp_global,    ONLY : intra_pool_comm
  USE constants,    ONLY : tpi, fpi
  USE gsmooth,      ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs, nls, ngms
  USE gvect,        ONLY : g, gg, gstart
  USE parameters,   ONLY : ntypx
  USE ions_base,    ONLY : nat, tau, atm, ityp

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  complex(dp), intent(in) :: rho_s(nrxxs)
  real(dp), intent(out) :: hfi_bare(nat), hfi_bare_zora(nat)  
  !-- local variables ----------------------------------------------------
#ifdef ZORA
  real(dp) ::  delta_Th(ngm, ntypx)
#endif
  integer :: ig, na
  real(dp) :: arg
  complex(dp) :: phase

  ! transform to reciprocal space
  call cft3s(rho_s, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)

#ifdef ZORA
  ! Fourier transform of Thomson's delta function
  call delta_thomson_radial_ft(delta_Th)
#endif
  
  ! fourier transform on the atomic position
  hfi_bare(:) = 0.d0
  hfi_bare_zora(:) = 0.d0
  do na = 1, nat
      do ig = gstart, ngms
          arg = sum(tau(1:3,na) * g(1:3,ig)) * tpi
          phase = cmplx(cos(arg),sin(arg), kind=dp)
          hfi_bare(na) = hfi_bare(na) + real(rho_s(nls(ig)) * phase, kind=dp)
#ifdef ZORA
          hfi_bare_zora(na) = hfi_bare_zora(na) + &
              real(delta_Th(ig,ityp(na)) * rho_s(nls(ig)) * phase, kind=dp)
#endif
      enddo
  enddo
  call mp_sum(hfi_bare, intra_pool_comm)
  call mp_sum(hfi_bare_zora, intra_pool_comm)

  return
END SUBROUTINE hfi_fc_bare_el



!-----------------------------------------------------------------------
SUBROUTINE delta_thomson_radial_ft(delta_th)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the Fourier transform of the Thomson-delta (for ZORA)
  !  
  USE kinds,                 ONLY : dp
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : ngm, gg
  USE ions_base,             ONLY : ntyp => nsp, atm
  USE constants,             ONLY : pi, fpi
  USE cell_base,             ONLY : tpiba
  USE gipaw_module,          ONLY : alpha, iverbosity, spline_integration, &
                                    radial_integral_splines
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: delta_th(ngm, ntyp)

  !-- local variables ----------------------------------------------------
  integer :: gv, j, nt
  real(dp) :: gr, r_Thomson
  real(dp), allocatable :: f_radial(:), work(:)
  integer, external :: atomic_number
  
  do nt = 1, ntyp
      allocate( work(rgrid(nt)%mesh), f_radial(rgrid(nt)%mesh) )
     
      ! Thomson's delta function
      r_Thomson = atomic_number(atm(nt)) * alpha ** 2
     
      ! terms rgrid(nt)%r(j) ** 2 from the definition of delta_Thomson
      ! and the radial volume element r^2 in integral cancel each other
      do j = 1, rgrid(nt)%mesh
          f_radial(j) = 2 / (fpi*r_Thomson * (1 + 2*rgrid(nt)%r(j)/r_Thomson)**2)
      enddo
     
      do gv = 1, ngm
          work = 0.0_dp
          do j = 1, rgrid(nt)%mesh
              gr = sqrt(gg(gv)) * tpiba * rgrid(nt)%r(j)
              if ( gr < 1.0d-8 ) then
                  work(j) = f_radial(j)*fpi
              else
                  work(j) = f_radial(j)*fpi * sin(gr)/gr
              endif
          enddo
          if (radial_integral_splines) then
              delta_th(gv,nt) = spline_integration(rgrid(nt)%r(:rgrid(nt)%mesh), &
                  work(:rgrid(nt)%mesh) )
          else
              call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab(:), delta_th(gv,nt))
          endif
     enddo

     deallocate (work, f_radial)
  enddo
  
END SUBROUTINE delta_thomson_radial_ft



!-----------------------------------------------------------------------
SUBROUTINE hfi_fc_gipaw_correction(fc_gipaw, fc_gipaw_zora)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the GIPAW contribution to the Fermi-contact
  !  
  USE kinds,                 ONLY : dp
  USE uspp,                  ONLY : ap
  USE parameters,            ONLY : lmaxx, ntypx
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : g,ngm,ecutwfc, gg
  USE klist,                 ONLY : nks, xk, wk
  USE cell_base,             ONLY : tpiba2
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp, atm
  USE wvfct,                 ONLY : npwx, nbnd, npw, igk, g2kin
  USE wavefunctions_module,  ONLY : evc
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE becmod,                ONLY : calbec
  USE constants,             ONLY : pi, fpi
  USE mp_global,             ONLY : intra_pool_comm
  USE mp,                    ONLY : mp_sum
  USE buffers,               ONLY : get_buffer
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE scf,                   ONLY : rho
  USE lsda_mod,              ONLY : current_spin, nspin, isk
  USE wvfct,                 ONLY : current_k, wg
  USE io_global,             ONLY : stdout
  USE gipaw_module,          ONLY : job, nbnd_occ, alpha, iverbosity, &
                                    spline_integration, &
                                    spline_integration_mirror, &
                                    radial_integral_splines
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: fc_gipaw(nat), fc_gipaw_zora(nat)
  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: at_hfi(:,:,:), at_hfi_zora(:,:,:)
  real(dp), allocatable :: work(:)
  integer :: j, nt, ibnd, il1, il2, ik, nbs1, nbs2, kkpsi
  integer :: m1, m2, lm1, lm2, l1, l2, nrc
  integer :: ijkb0, ih, jh, na, ikb, jkb
  integer :: s_min, s_maj, s_weight, r_first
  real(dp) :: r_Thomson
  complex(dp) :: bec_product
  integer, external :: atomic_number
  
  allocate( at_hfi(paw_nkb,paw_nkb,ntyp) )
  allocate( at_hfi_zora(paw_nkb,paw_nkb,ntyp) )
  at_hfi = 0.0_dp
  at_hfi_zora = 0.0_dp
  
  ! calculate radial integrals: <aephi|aephi> - <psphi|psphi>
  do nt = 1, ntyp
     kkpsi = paw_recon(nt)%aephi(1)%kkpsi
     allocate( work(kkpsi) )

     r_first = 1
     if ( abs ( rgrid(nt)%r(1) ) < 1d-8 ) r_first = 2
     r_thomson = atomic_number(atm(nt)) * alpha**2
     
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
           at_hfi(il1,il2,nt) = work(r_first)
           
#ifdef ZORA           
           ! multiply with Thomson's delta function
           do j = r_first, nrc
              work(j) = work(j) * 2/(r_thomson * (1+2*rgrid(nt)%r(j)/r_thomson)**2)
           enddo
           
           if (radial_integral_splines) then
              at_hfi_zora(il1,il2,nt) = spline_integration(rgrid(nt)%r(:nrc), work(:nrc))
           else
              call simpson(nrc, work, rgrid(nt)%rab(:), at_hfi_zora(il1,il2,nt))
           endif
#endif
        enddo
     enddo

     deallocate ( work )
  enddo
  
  !  calculate the reconstruction part
  fc_gipaw = 0.d0
  fc_gipaw_zora = 0.d0
  call select_spin(s_min, s_maj)  
  
  do ik = 1, nks
     current_k = ik
     current_spin = isk(ik)
     
     if (current_spin == s_min) then
        s_weight = -1
     else
        s_weight = +1
     endif
     
     call gk_sort ( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
     call get_buffer ( evc, nwordwfc, iunwfc, ik)
     
     call init_gipaw_2 ( npw, igk, xk(1,ik), paw_vkb )
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
                       
                       fc_gipaw(na) = fc_gipaw(na) + s_weight * at_hfi(nbs1,nbs2,nt) &
                            * bec_product * wg(ibnd,ik)
#ifdef ZORA                       
                       fc_gipaw_zora(na) = fc_gipaw_zora(na) + s_weight * at_hfi_zora(nbs1,nbs2,nt) &
                            * bec_product * wg(ibnd,ik)
#endif                       
                    enddo
                 enddo
                 
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              endif
              
           enddo
        enddo
     enddo
     
  enddo
  
END SUBROUTINE hfi_fc_gipaw_correction


