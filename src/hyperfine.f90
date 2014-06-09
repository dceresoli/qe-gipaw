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
  USE kinds,                  ONLY : dp 
  USE io_global,              ONLY : stdout
  USE parameters,             ONLY : ntypx
  USE constants,              ONLY : pi, tpi, fpi, angstrom_au, rytoev, electronvolt_si, c_si
  USE fft_base,               ONLY : dffts, dfftp
  USE scf,                    ONLY : rho
  USE symme,                  ONLY : symtensor
  USE lsda_mod,               ONLY : current_spin, nspin
  USE wvfct,                  ONLY : current_k
  USE ions_base,              ONLY : nat, tau, atm, ityp
  use constants,              ONLY : bohr_radius_si
  USE mp_pools,               ONLY : intra_pool_comm
  USE mp,                     ONLY : mp_sum
  USE gipaw_module,           ONLY : hfi_nuclear_g_factor, hfi_output_unit, &
                                     job, iverbosity, core_relax_method, &
                                     hfi_via_reconstruction_only
 
  !-- constants ----------------------------------------------------------
  IMPLICIT NONE
  real(dp), parameter :: mu0_by_fpi = 1e-7
  real(dp), parameter :: mu_n = 5.05078324e-27_dp
  real(dp), parameter :: bohr_radius = bohr_radius_si
  real(dp), parameter :: gamma_e = 28024953.64_dp
  real(dp), parameter :: lambda = C_SI / 1.0e+8_dp
  real(dp), parameter :: common_factor = mu0_by_fpi * mu_n / Bohr_radius ** 3

  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: spin_den(:), rho_s(:,:)
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
  allocate( rho_s(dffts%nnr,2), spin_den(dffts%nnr) )
  call get_smooth_density(rho_s)  ! this subroutine is in efg.g90
  spin_den(:) = rho_s(:,s_maj) - rho_s(:,s_min)
  call efg_bare_el(spin_den, hfi_dip_bare)

  ! calculate GIPAW dipole correction
  call efg_correction(hfi_dip_gipaw)


  !--------------------------------------------------------------------
  ! Fermi-contact contribution to hyperfine
  !--------------------------------------------------------------------
  allocate( hfi_fc_bare(nat), hfi_fc_bare_zora(nat) )
  allocate( hfi_fc_gipaw(nat), hfi_fc_gipaw_zora(nat) )
  allocate( hfi_fc_core(nat), hfi_fc_tot(nat) )

  ! calculate the bare Fermi-contact contribution
  spin_den(:) = rho_s(:,s_maj) - rho_s(:,s_min)
  call hfi_fc_bare_el(spin_den, hfi_fc_bare, hfi_fc_bare_zora)
  deallocate( rho_s, spin_den )

  ! calculate the GIPAW Fermi-contact contribution
  call hfi_fc_gipaw_correction(hfi_fc_gipaw, hfi_fc_gipaw_zora)

  ! calculate the core-relaxation Fermi-contact contribution
  call hfi_fc_core_relax(core_relax_method, hfi_fc_core)


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

  ! Print the Fermi contact term
  write(stdout,*)
  write(stdout,'(5X,''ISOTROPIC (FERMI-CONTACT) COUPLINGS WITHOUT ZORA:'')')
  write(stdout,'(5X,''USING CORE-RELAXATION METHOD: PRB 76, 035124 (2007)'')')
  write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
  if (iverbosity > 1) then
    write(stdout,'(5X,''----- spin-densities in bohrradius^-3 -----'')')
    write(stdout,'(5X,8X,''  bare            GIPAW           core-relax      total'')')
    do na = 1, nat
        hfi_fc_tot(na) = hfi_fc_bare(na) + hfi_fc_gipaw(na) + hfi_fc_core(na)
        write(stdout,1002) atm(ityp(na)), na, hfi_fc_bare(na), &
            hfi_fc_gipaw(na), hfi_fc_core(na), hfi_fc_tot(na)
    enddo
  endif
  write(stdout,'(5X,''----- Fermi contact in '',A,'' -----'')') trim(hfi_output_unit)
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
!  write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
  write(stdout,'(5X,''Warning: core-relaxation is left out here. For the corrections see above.'')')
  write(stdout,'(5X,8X,''  bare            GIPAW           core-relax      total'')')
  do na = 1, nat
      hfi_fc_tot(na) = hfi_fc_bare_zora(na) + hfi_fc_gipaw_zora(na) ! + hfi_fc_core(na)
      fact = 8*pi/3 * output_factor * hfi_nuclear_g_factor(ityp(na))
      write(stdout,1002) atm(ityp(na)), na, fact * hfi_fc_bare_zora(na), &
!          fact * hfi_fc_gipaw_zora(na), fact * hfi_fc_core(na), fact * hfi_fc_tot(na)
          fact * hfi_fc_gipaw_zora(na), fact * 0.d0, fact * hfi_fc_tot(na)
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
  USE kinds,                  ONLY : dp 
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : intra_pool_comm
  USE constants,              ONLY : tpi, fpi
  USE gvecs,                  ONLY : nls, ngms
  USE gvect,                  ONLY : g, gg, gstart
  USE parameters,             ONLY : ntypx
  USE ions_base,              ONLY : nat, tau, atm, ityp
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : fwfft
  USE gipaw_module,           ONLY : hfi_via_reconstruction_only

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(in) :: rho_s(dffts%nnr)
  real(dp), intent(out) :: hfi_bare(nat), hfi_bare_zora(nat)  
  !-- local variables ----------------------------------------------------
#ifdef ZORA
  real(dp) ::  delta_Th(ngm, ntypx)
#endif
  complex(dp), allocatable :: rhoaux(:)
  integer :: ig, na
  real(dp) :: arg
  complex(dp) :: phase

  hfi_bare(:) = 0.d0
  hfi_bare_zora(:) = 0.d0
  IF ( hfi_via_reconstruction_only ) RETURN

  ! transform to reciprocal space
  allocate(rhoaux(dffts%nnr))
  rhoaux(:) = cmplx(rho_s(:), kind=dp)
  CALL fwfft('Smooth', rhoaux, dffts)

#ifdef ZORA
  ! Fourier transform of Thomson's delta function
  call delta_thomson_radial_ft(delta_Th)
#endif
  
  ! fourier transform on the atomic position
  do na = 1, nat
      do ig = gstart, ngms
          arg = sum(tau(1:3,na) * g(1:3,ig)) * tpi
          phase = cmplx(cos(arg),sin(arg), kind=dp)
          hfi_bare(na) = hfi_bare(na) + real(rhoaux(nls(ig)) * phase, kind=dp)
#ifdef ZORA
          hfi_bare_zora(na) = hfi_bare_zora(na) + &
              real(delta_Th(ig,ityp(na)) * rhoaux(nls(ig)) * phase, kind=dp)
#endif
      enddo
  enddo
  call mp_sum(hfi_bare, intra_pool_comm)
  call mp_sum(hfi_bare_zora, intra_pool_comm)

  deallocate(rhoaux)
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
  USE gipaw_module,          ONLY : alpha, iverbosity
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: delta_th(ngm, ntyp)

  !-- local variables ----------------------------------------------------
  integer :: gv, j, nt
  real(dp) :: gr, r_Thomson, norm_Th(ngm,ntyp) 
  real(dp), allocatable :: f_radial(:), work(:), work_norm(:)
  integer, external :: atomic_number
  
  do nt = 1, ntyp
      allocate( work(rgrid(nt)%mesh), f_radial(rgrid(nt)%mesh) )
      allocate( work_norm(rgrid(nt)%mesh) )
     
      ! Thomson's delta function
      r_Thomson = atomic_number(atm(nt)) * alpha ** 2
     
      ! terms rgrid(nt)%r(j) ** 2 from the definition of delta_Thomson
      ! and the radial volume element r^2 in integral cancel each other
      do j = 1, rgrid(nt)%mesh
          f_radial(j) = 2 / (fpi*r_Thomson * (1 + 2*rgrid(nt)%r(j)/r_Thomson)**2)
      enddo
     
      do gv = 1, ngm
          work = 0.0_dp
          work_norm = 0.0_dp
          do j = 1, rgrid(nt)%mesh
              work_norm(j) = f_radial(j) * fpi
              gr = sqrt(gg(gv)) * tpiba * rgrid(nt)%r(j)
              if ( gr < 1.0d-8 ) then
                  work(j) = f_radial(j)*fpi
              else
                  work(j) = f_radial(j)*fpi * sin(gr)/gr
              endif
          enddo
          call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab(:), delta_th(gv,nt))
          call simpson(rgrid(nt)%mesh, work_norm, rgrid(nt)%rab(:), norm_th(gv,nt))
     enddo
!    correction od the norm 

     DO gv = 1, ngm
        delta_th(gv,nt) = delta_th(gv,nt)/norm_th(gv,nt)
     END DO
     deallocate (work_norm)

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
  USE gvect,                 ONLY : g, ngm, gg
  USE klist,                 ONLY : nks, xk, wk
  USE cell_base,             ONLY : tpiba2
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp, atm
  USE wvfct,                 ONLY : npwx, nbnd, npw, igk, g2kin, ecutwfc
  USE wavefunctions_module,  ONLY : evc
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE becmod,                ONLY : calbec
  USE constants,             ONLY : pi, fpi
  USE mp_pools,              ONLY : intra_pool_comm, inter_pool_comm
  USE mp,                    ONLY : mp_sum
  USE buffers,               ONLY : get_buffer
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE scf,                   ONLY : rho
  USE lsda_mod,              ONLY : current_spin, nspin, isk
  USE wvfct,                 ONLY : current_k, wg
  USE io_global,             ONLY : stdout
  USE gipaw_module,          ONLY : job, nbnd_occ, alpha, iverbosity, hfi_via_reconstruction_only

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: fc_gipaw(nat), fc_gipaw_zora(nat)
  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: at_hfi(:,:,:), at_hfi_zora(:,:,:), at_hfi_extra(:,:,:)
  real(dp), allocatable :: work(:)
!
  REAL ( dp ), ALLOCATABLE :: x_extrapolate(:), y_extrapolate(:)
  REAL ( dp ) :: x
  INTEGER :: hfi_extrapolation_npoints = 10000
  INTEGER :: norder_extrapolate = 3
!
  integer :: j, nt, ibnd, il1, il2, ik, nbs1, nbs2, kkpsi
  integer :: m1, m2, lm1, lm2, l1, l2, nrc
  integer :: ijkb0, ih, jh, na, ikb, jkb
  integer :: s_min, s_maj, s_weight, r_first
  real(dp) :: r_Thomson
  complex(dp) :: bec_product
  integer, external :: atomic_number
  
  allocate( at_hfi(paw_nkb,paw_nkb,ntyp) )
  allocate( at_hfi_zora(paw_nkb,paw_nkb,ntyp) )
  allocate( at_hfi_extra(paw_nkb,paw_nkb,ntyp) )
  at_hfi = 0.0_dp
  at_hfi_zora = 0.0_dp
  at_hfi_extra = 0.0_dp

  ALLOCATE ( x_extrapolate(hfi_extrapolation_npoints) )
  ALLOCATE ( y_extrapolate(hfi_extrapolation_npoints) )
  
  ! calculate radial integrals: <aephi|aephi> - <psphi|psphi>
  do nt = 1, ntyp
     kkpsi = paw_recon(nt)%aephi(1)%kkpsi
     allocate( work(kkpsi) )

     r_first = 1
     if ( abs ( rgrid(nt)%r(1) ) < 1d-8 ) r_first = 2
     r_thomson = atomic_number(atm(nt)) * alpha**2
     


     DO j = 1, hfi_extrapolation_npoints
        x_extrapolate(j) = j / REAL ( hfi_extrapolation_npoints + 1, dp ) &
             * rgrid(nt)%r(r_first)
     END DO



     do il1 = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(il1)%label%nrc
        l1 = paw_recon(nt)%psphi(il1)%label%l
        if (l1 /= 0) cycle
        
        do il2 = 1, paw_recon(nt)%paw_nbeta
           l2 = paw_recon(nt)%psphi(il2)%label%l
           if (l2 /= 0) cycle
           
           work = 0.0_dp
           IF ( hfi_via_reconstruction_only ) THEN
              do j = r_first, nrc
                 work(j) = &
                      ( paw_recon(nt)%aephi(il1)%psi(j) &
                      * paw_recon(nt)%aephi(il2)%psi(j) ) &
                      / rgrid(nt)%r(j) ** 2 / fpi
              end do
           ELSE
              do j = r_first, nrc
                 work(j) = &
                      ( paw_recon(nt)%aephi(il1)%psi(j) &
                      * paw_recon(nt)%aephi(il2)%psi(j) &
                      - paw_recon(nt)%psphi(il1)%psi(j) &
                      * paw_recon(nt)%psphi(il2)%psi(j) ) &
                      / rgrid(nt)%r(j) ** 2 / fpi
              enddo
           END IF
           
           ! density at the origin
           at_hfi(il1,il2,nt) = work(r_first)
           
#ifdef ZORA
           ! For ZORA (pseudos in scalar-relativistic approximation) we need to extrapolate:
           x=0.d0
           at_hfi_extra(il1,il2,nt) &
                = spline_mirror_extrapolate ( rgrid(nt)%r(:nrc), work(:nrc), x )

           CALL radial_extrapolation ( rgrid(nt)%r(:), work(:nrc), &
                x_extrapolate, y_extrapolate, norder_extrapolate )

           
           ! multiply with Thomson's delta function
           do j = r_first, nrc
              work(j) = work(j) * 2/(r_thomson * (1+2*rgrid(nt)%r(j)/r_thomson)**2)
           enddo
           
           call simpson(nrc, work, rgrid(nt)%rab(:), at_hfi_zora(il1,il2,nt))

           do j = 1, hfi_extrapolation_npoints
              y_extrapolate(j) = y_extrapolate(j) &
                   * 2 / ( r_thomson &
                   * ( 1 + 2 * x_extrapolate(j) / r_thomson ) ** 2 )
           end do

           at_hfi_extra(il1,il2,nt) = at_hfi_zora(il1,il2,nt) &
                   + spline_integration_mirror ( x_extrapolate, y_extrapolate )
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
!                       fc_gipaw_zora(na) = fc_gipaw_zora(na) + s_weight * at_hfi_zora(nbs1,nbs2,nt) &
                       fc_gipaw_zora(na) = fc_gipaw_zora(na) + s_weight * at_hfi_extra(nbs1,nbs2,nt) &
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

  call mp_sum( fc_gipaw, inter_pool_comm )
  call mp_sum( fc_gipaw_zora, inter_pool_comm )
  
  deallocate( at_hfi )
  deallocate( at_hfi_zora )
  deallocate( at_hfi_extra )




CONTAINS

  FUNCTION spline_integration_mirror ( xdata, ydata )

    !                                                                                                                                                 
    ! Like 'spline_integration' but assumes the function to be symmetric                                                                              
    !    on x axis [i.e. f(-x) = f(x) ]                                                                                                               

    USE splinelib, ONLY : spline

    IMPLICIT NONE

    ! Return type                                                                                                                                     
    REAL ( dp ) :: spline_integration_mirror

    ! Arguments                                                                                                                                       
    REAL ( dp ), INTENT ( IN ) :: xdata(:), ydata(:)

    ! Local                                                                                                                                           
    INTEGER  :: n
    REAL ( dp ) :: startd, startu
    REAL ( dp ) :: xdata2(2*SIZE(xdata)), ydata2(2*SIZE(ydata))
    REAL ( dp ) :: d2y(2*SIZE(xdata))

    !--------------------------------------------------------------------------                                                                       

    IF ( SIZE ( xdata ) /= SIZE ( ydata ) ) &
         CALL errore ( "spline_interpolation", &
         "sizes of arguments do not match", 1 )

    xdata2(SIZE ( xdata ):1:-1) = - xdata
    xdata2(SIZE ( xdata )+1:) = xdata

    ydata2(SIZE ( xdata ):1:-1) = ydata
    ydata2(SIZE ( xdata )+1:) = ydata

    startu = ( ydata(2) - ydata(1) ) / ( xdata(2) - xdata(1) )
    startu = 0.0
    startd = 0.0

    CALL spline ( xdata2, ydata2, startu, startd, d2y )

    spline_integration_mirror = 0.0
    DO n = 1, SIZE ( xdata2 )
       spline_integration_mirror = spline_integration_mirror &
            + splint_integr ( xdata2, ydata2, d2y, xdata2(n) )
    END DO

    spline_integration_mirror = spline_integration_mirror / 2

  END FUNCTION spline_integration_mirror





  FUNCTION spline_mirror_extrapolate ( xdata, ydata, x )

    USE splinelib, ONLY : spline, splint

    IMPLICIT NONE

    ! Return type                                                                                                                                     
    REAL ( dp ) :: spline_mirror_extrapolate

    ! Arguments                                                                                                                                       
    REAL ( dp ), INTENT ( IN ) :: xdata(:), ydata(:), x

    ! Local                                                                                                                                           
    INTEGER  :: n
    REAL ( dp ) :: startd, startu
    REAL ( dp ) :: xdata2(2*SIZE(xdata)), ydata2(2*SIZE(ydata))
    REAL ( dp ) :: d2y(2*SIZE(xdata))

    !--------------------------------------------------------------------------                                   

    IF ( SIZE ( xdata ) /= SIZE ( ydata ) ) &
         CALL errore ( "spline_mirror_extrapolate", &
         "sizes of arguments do not match", 1 )

    n = SIZE ( xdata )

    xdata2(n:1:-1) = - xdata
    xdata2(n+1:) = xdata

    ydata2(n:1:-1) = ydata
    ydata2(n+1:) = ydata

    startu = 0.0
    startd = 0.0

    CALL spline ( xdata2, ydata2, startu, startd, d2y )

    spline_mirror_extrapolate = splint ( xdata2, ydata2, d2y, x )

  END FUNCTION spline_mirror_extrapolate



  !****************************************************************************                                                                        

  FUNCTION splint_integr ( xdata, ydata, d2y, x )

    USE splinelib, ONLY : spline

    IMPLICIT NONE

    ! Return type                                                                                                                                      
    REAL ( dp ) :: splint_integr

    ! Arguments                                                                                                                                        
    REAL ( dp ), INTENT ( IN ) :: xdata(:), ydata(:), d2y(:), x

    ! Local                                                                                                                                            
    INTEGER  :: khi, klo, xdim
    REAL ( dp ) :: a, b, da, db, h

    !--------------------------------------------------------------------------                                                                        

    xdim = SIZE( xdata )

    klo = 1
    khi = xdim

    klo = MAX( MIN( locate( xdata, x ), ( xdim - 1 ) ), 1 )

    khi = klo + 1

    h = xdata(khi) - xdata(klo)

    a = ( xdata(khi) - x ) / h
    b = ( x - xdata(klo) ) / h
    da = -1 / h
    db =  1 / h

    splint_integr = -0.5 * ydata(klo) / da + 0.5 * ydata(khi) / db &
         + (  0.25 / da * d2y(klo) &
         + ( -0.25 / db * d2y(khi) ) ) &
         * ( h**2 ) / 6.D0

  END FUNCTION splint_integr



         !-------------------------------------------------------------------                                                                          
         FUNCTION locate( xx, x )
           !-------------------------------------------------------------------                                                                        
           !                                                                                                                                           
           IMPLICIT NONE
           !                                                                                                                                           
           REAL(DP), INTENT(IN) :: xx(:)
           REAL(DP), INTENT(IN) :: x
           !                                                                                                                                           
           INTEGER :: locate
           INTEGER :: n, jl, jm, ju
           LOGICAL :: ascnd
           !                                                                                                                                           
           !                                                                                                                                           
           n     = SIZE( xx )
           ascnd = ( xx(n) >= xx(1) )
           jl    = 0
           ju    = n + 1
           !                                                                                                                                           
           main_loop: DO
              !                                                                                                                                        
              IF ( ( ju - jl ) <= 1 ) EXIT main_loop
              !                                                                                                                                        
              jm = ( ju + jl ) / 2
              !                                                                                                                                        
              IF ( ascnd .EQV. ( x >= xx(jm) ) ) THEN
                 !                                                                                                                                     
                 jl = jm
                 !                                                                                                                                     
              ELSE
                 !                                                                                                                                     
                 ju = jm
                 !                                                                                                                                     
              END IF
              !                                                                                                                                        
           END DO main_loop
           !                                                                                                                                           
           IF ( x == xx(1) ) THEN
              !                                                                                                                                        
              locate = 1
              !                                                                                                                                        
           ELSE IF ( x == xx(n) ) THEN
              !                                                                                                                                        
              locate = n - 1
              !                                                                                                                                        
           ELSE
              !                                                                                                                                        
              locate = jl
              !                                                                                                                                        
           END IF
           !                                                                                                                                           
         END FUNCTION locate

  !****************************************************************************                                                                        

  SUBROUTINE radial_extrapolation ( x, y, x_extrapolate, y_extrapolate, &
       norders )

    IMPLICIT NONE

    ! Arguments                                                                                                                                        
    REAL ( dp ), INTENT ( IN ) :: x(:), y(:)
    REAL ( dp ), INTENT ( IN ) :: x_extrapolate(:)
    REAL ( dp ), INTENT ( OUT ) :: y_extrapolate(:)
    INTEGER, INTENT ( IN ) :: norders

    ! Local                                                                                                                                            
    INTEGER :: n, i, j

    REAL ( dp ) :: a(0:norders,0:norders), b(0:norders), c(0:norders)

    !--------------------------------------------------------------------------                                                                        

    ! Dirac delta function                                                                                                                             
    do j = 0, norders
       a(0:norders,j) = x(1:norders+1) ** j
       b(0:norders) = y(1:norders+1)
    END DO

    CALL invert ( a, norders + 1 )

    c = MATMUL ( a, b )

    y_extrapolate = 0.0
    DO i = 1, SIZE ( x_extrapolate )
       DO j = 0, norders
          y_extrapolate(i) = y_extrapolate(i) + c(j) * x_extrapolate(i) ** j
       END DO
    END DO

  END SUBROUTINE radial_extrapolation

  !****************************************************************************                                         


  SUBROUTINE invert ( a, n )

    IMPLICIT NONE

    ! Arguments                                                                                                                                        
    INTEGER, INTENT ( IN ) :: n
    REAL ( dp ), INTENT ( INOUT ) :: a(n,n)

    ! Local                                                                                                                                            
    INTEGER :: ipvt(n)
    INTEGER :: info
    REAL ( dp ) :: cwork(n)

    !--------------------------------------------------------------------------                                                                        

    CALL dgetrf ( n, n, a, n, ipvt, info )
    IF ( info /= 0 ) &
         CALL errore ( "invert_in_hypefile", "dgetrf failed", ABS ( info ) )

    CALL dgetri ( n, a, n, ipvt, cwork, n, info )
    IF ( info /= 0 ) &
         CALL errore ( "invert_in_hypefile", "dgetri failed", ABS ( info ) )

  END SUBROUTINE invert

  
END SUBROUTINE hfi_fc_gipaw_correction





