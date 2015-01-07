!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE gipaw_setup
  !-----------------------------------------------------------------------
  !
  ! ... GIPAW setup
  !
  USE kinds,         ONLY : dp
  USE io_global,     ONLY : stdout, ionode
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, atm
  USE atom,          ONLY : rgrid
  USE wvfct,         ONLY : nbnd, et, wg, npwx
  USE lsda_mod,      ONLY : nspin, lsda
  USE scf,           ONLY : v, vrs, vltot, rho, rho_core, kedtau
  USE gvect,         ONLY : ngm
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE klist,         ONLY : xk, degauss, ngauss, nks, nelec, lgauss, wk, two_fermi_energies
  USE ktetra,        ONLY : ltetra
  USE noncollin_module,  ONLY : noncolin
  USE constants,     ONLY : degspin, pi
  USE symm_base,     ONLY : nsym, s
  USE mp_pools,      ONLY : inter_pool_comm 
  USE mp,            ONLY : mp_max, mp_min 
  USE dfunct,        ONLY : newd
  USE pwcom,         ONLY : ef
  USE constants,     ONLY : rytoev
  USE gipaw_module

  implicit none
  integer :: ik, ibnd, ipol
  real(dp) :: emin, emax, xmax, small, fac, target
    
  call start_clock ('gipaw_setup')
    
  ! TODO: test whether the symmetry operations map the Cartesian axis to each
  ! call test_symmetries ( s, nsym )    

  ! initialize pseudopotentials and projectors for LDA+U
  call init_us_1
  call init_at_1

  ! setup GIPAW operators
  call gipaw_setup_integrals
  call gipaw_setup_l

  ! computes the total local potential (external+scf) on the smooth grid
  call setlocal
  call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
    
  ! compute the D for the pseudopotentials
  call newd
    
  !! set non linear core correction stuff (IS THIS REALLY NEEDED?)
  !! nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !!if (nlcc_any) allocate (drc( ngm, ntyp))
  !! setup all gradient correction stuff
  !!call setup_dgc

  ! some pre-conditions
  if (ltetra) call errore('gipaw_setup','GIPAW + tetrahedra not implemented', 1)
  if (noncolin) call errore('gipaw_setup','GIPAW + non-collinear not implemented', 1)
  if (two_fermi_energies .and. lgauss) &
     call errore('gipaw_setup','GIPAW + two Fermi energies not implemented', 1)

  ! computes the number of occupied bands for each k point
  nbnd_occ (:) = 0
  if (lgauss) then
     write(stdout,*)
     write(stdout,'(5X,''smearing ngauss='',I4,2X,''degauss='',F8.4,'' Ry'')') &
          ngauss, degauss
     ! discard conduction bands such that w0gauss(x,n) < small
     ! hint:
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     small = 6.3491173359333d-8

     ! appropriate limit for gaussian broadening (used for all ngauss)
     xmax = sqrt(-log(sqrt(pi)*small))

     ! appropriate limit for Fermi-Dirac
     if (ngauss == -99) then
        fac = 1.d0 / sqrt(small)
        xmax = 2.d0 * log(0.5d0*(fac + sqrt(fac*fac-4.d0)))
     endif
     target = ef + xmax * degauss
     do ik = 1, nks
        do ibnd = 1, nbnd
!DEBUG           if (ionode) write(70,*) et(ibnd,ik), wg(ibnd,ik)/wk(ik)
           if (et(ibnd,ik) < target) nbnd_occ(ik) = ibnd
        enddo
        if (nbnd_occ (ik) == nbnd) &
           write(stdout,'(5X,''Possibly too few bands at k-point:'',I6)') ik
     enddo
  else 
    ! general case
     do ik = 1, nks
       do ibnd = 1, nbnd
         if (wk(ik) > 0.d0) then
           if (wg(ibnd,ik)/wk(ik) > 1d-4 ) nbnd_occ(ik) = ibnd
          endif
       end do
     end do
  end if
    
  ! computes alpha_pv
  emin = et (1, 1)
  do ik = 1, nks
    do ibnd = 1, nbnd
      emin = min (emin, et (ibnd, ik) )
    enddo
  enddo
#ifdef __MPI
  ! find the minimum across pools
  call mp_min( emin, inter_pool_comm )
#endif

  if (lgauss) then
     ! metal
     emax = target
     alpha_pv = emax - emin
  else
     ! insulator
     emax = et(1,1)
     do ik = 1, nks
        do ibnd = 1, nbnd_occ(ik)
           emax = max(emax, et(ibnd,ik))
        enddo
     enddo
#ifdef __MPI
     ! find the maximum across pools
     call mp_max( emax, inter_pool_comm )
#endif
     alpha_pv = 2.d0 * (emax - emin)
  endif

  ! avoid zero value for alpha_pv
  alpha_pv = max(alpha_pv, 1.0d-2)
  write(stdout,'(5X,''alpha_pv='',F12.4,'' eV'')') alpha_pv*rytoev

  if (iverbosity > 10) then
     write(stdout,*)
     write(stdout,'(5X,''Number of occupied bands for each k-point:'')')
     do ik = 1, nks
        write(stdout,'(5X,''k-point:'',I6,4X,''nbnd_occ='',I4)') ik, nbnd_occ(ik)
     enddo
     write(stdout,*)
  endif

  call stop_clock('gipaw_setup')
    
END SUBROUTINE gipaw_setup



!-----------------------------------------------------------------------
SUBROUTINE gipaw_setup_integrals
  !-----------------------------------------------------------------------
  !
  ! ... Setup the GIPAW integrals: NMR core contribution, diamagnetic 'E'
  ! ... and paramagnetic 'F' terms, relativistic mass corrections
  !
  USE gipaw_module
  USE kinds,         ONLY : dp
  USE io_global,     ONLY : stdout, ionode
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, atm
  USE atom,          ONLY : rgrid
  USE paw_gipaw,     ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp, set_paw_upf
  USE uspp_param,    ONLY : upf
  USE mp_pools,      ONLY : inter_pool_comm 
  USE wvfct,         ONLY : nbnd, npwx

  implicit none

  real(dp), allocatable :: work(:), kinetic_aephi(:), kinetic_psphi(:)
  real(dp), allocatable :: aephi_dvloc_dr(:), psphi_dvloc_dr(:)
  integer :: nt, il, il1, il2, l1, l2, j, kkpsi, nrc
  integer :: core_orb
  real(dp) :: integral, occupation
    
  ! initialize data, also for the case that no GIPAW is present
  if ( .not. allocated(paw_recon) ) allocate(paw_recon(ntyp))
  paw_recon(:)%gipaw_data_in_upf_file = .false.
  paw_recon(:)%paw_nbeta = 0
  paw_recon(:)%paw_nh = 0
  paw_recon(:)%paw_kkbeta = 0
  paw_recon(:)%gipaw_ncore_orbital = 0
  paw_recon(:)%vloc_present = .false.
  do nt = 1, ntyp
    paw_recon(nt)%paw_lll(:) = 0
  end do
    

  ! setup GIPAW projectors
  do nt = 1, ntyp
    call set_paw_upf(nt, upf(nt))
    !!call read_recon(file_reconstruction(nt), nt, paw_recon(nt))
  enddo
    
  do nt = 1, ntyp
    do il = 1, paw_recon(nt)%paw_nbeta
      if ( paw_recon(nt)%psphi(il)%label%rc < -0.99d0 ) then
        rc(nt,paw_recon(nt)%psphi(il)%label%l) = 1.6d0
        rc(nt,paw_recon(nt)%aephi(il)%label%l) = 1.6d0
        paw_recon(nt)%psphi(il)%label%rc = rc(nt,paw_recon(nt)%psphi(il)%label%l)
        paw_recon(nt)%aephi(il)%label%rc = rc(nt,paw_recon(nt)%aephi(il)%label%l)
      else
        rc(nt,paw_recon(nt)%psphi(il)%label%l) = paw_recon(nt)%psphi(il)%label%rc
        rc(nt,paw_recon(nt)%aephi(il)%label%l) = paw_recon(nt)%aephi(il)%label%rc
      endif
    enddo
  enddo
    
  call init_gipaw_1()
 
  ! allocate GIPAW projectors   
  allocate ( paw_vkb(npwx,paw_nkb) )
  allocate ( paw_becp(paw_nkb,nbnd) )
  allocate ( paw_becp2(paw_nkb,nbnd) )
  allocate ( paw_becp3(paw_nkb,nbnd) )

  ! allocate GIPAW integrals    
  allocate ( radial_integral_diamagnetic(nbrx,nbrx,ntypx) )
  allocate ( radial_integral_paramagnetic(nbrx,nbrx,ntypx) )
  allocate ( radial_integral_diamagnetic_so(nbrx,nbrx,ntypx) )
  allocate ( radial_integral_paramagnetic_so(nbrx,nbrx,ntypx) )
  allocate ( radial_integral_rmc(nbrx,nbrx,ntypx) )
  radial_integral_diamagnetic = 0.d0
  radial_integral_paramagnetic = 0.d0
  radial_integral_diamagnetic_so = 0.d0
  radial_integral_paramagnetic_so = 0.d0
  radial_integral_rmc = 0.d0
    

  ! calculate GIPAW integrals
  do nt = 1, ntyp

    do il1 = 1, paw_recon(nt)%paw_nbeta
      l1 = paw_recon(nt)%psphi(il1)%label%l

      kkpsi = paw_recon(nt)%aephi(il1)%kkpsi
      nrc = paw_recon(nt)%psphi(il1)%label%nrc
          
      allocate ( work(kkpsi) )
          
      do il2 = 1, paw_recon(nt)%paw_nbeta
        l2 = paw_recon(nt)%psphi(il2)%label%l
            
        if ( l1 /= l2 ) cycle
             
        ! NMR diamagnetic: (1/r)
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / rgrid(nt)%r(j)
        enddo
        call simpson( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_diamagnetic(il1,il2,nt) )
             
        ! NMR paramagnetic: (1/r^3)
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     / rgrid(nt)%r(j) ** 3
        enddo
        call simpson( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_paramagnetic(il1,il2,nt) )
             
        ! calculate the radial integral only if the radial potential is present
        if ( .not. paw_recon(nt)%vloc_present ) cycle
             
        ! g-tensor relativistic mass correction: (-nabla^2)
        allocate ( kinetic_aephi ( kkpsi ), kinetic_psphi ( kkpsi ) )
        call radial_kinetic_energy (nrc, l2, rgrid(nt)%r(:nrc), paw_recon(nt)%aephi(il2)%psi(:nrc), &
                                    kinetic_aephi(:nrc))
        call radial_kinetic_energy (nrc, l2, rgrid(nt)%r(:nrc), paw_recon(nt)%psphi(il2)%psi(:nrc), &
                                    kinetic_psphi(:nrc))
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * kinetic_aephi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * kinetic_psphi(j) )
        enddo
        deallocate ( kinetic_aephi, kinetic_psphi )
        call simpson ( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_rmc(il1,il2,nt) )
             
        ! calculate dV/dr
        allocate ( aephi_dvloc_dr(nrc), psphi_dvloc_dr(nrc) )
        call radial_derivative (nrc, rgrid(nt)%r(:nrc), paw_recon(nt)%gipaw_ae_vloc(:nrc), aephi_dvloc_dr(:nrc))
        call radial_derivative (nrc, rgrid(nt)%r(:nrc), paw_recon(nt)%gipaw_ps_vloc(:nrc), psphi_dvloc_dr(:nrc))
             
        ! g tensor diamagnetic: (r*dV/dr)
        do j = 1, nrc
           work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                     - paw_recon(nt)%psphi(il1)%psi(j) * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                     * rgrid(nt)%r(j)
        enddo
        call simpson( nrc, work, rgrid(nt)%rab(:nrc), radial_integral_diamagnetic_so(il1,il2,nt) )
             
        ! g tensor paramagnetic: (1/r*dV/dr)
        do j = 1, nrc
          work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) * aephi_dvloc_dr(j) * paw_recon(nt)%aephi(il2)%psi(j) &
                    - paw_recon(nt)%psphi(il1)%psi(j) * psphi_dvloc_dr(j) * paw_recon(nt)%psphi(il2)%psi(j) ) &
                    / rgrid(nt)%r(j)
        enddo
        call simpson( nrc,work,rgrid(nt)%rab(:nrc), radial_integral_paramagnetic_so(il1,il2,nt) )
             
        deallocate ( aephi_dvloc_dr, psphi_dvloc_dr )
             
      enddo  ! l2
          
      deallocate ( work )
          
     enddo  ! l1
  enddo  ! nt


  ! Compute the shift due to core orbitals (purely diamagnetic)
  do nt = 1, ntyp
    if ( paw_recon(nt)%gipaw_ncore_orbital == 0 ) cycle

    allocate ( work(rgrid(nt)%mesh) )
    nmr_shift_core(nt) = 0.0
       
    do core_orb = 1, paw_recon(nt)%gipaw_ncore_orbital
      do j = 1, size(work)
        work(j) = paw_recon(nt)%gipaw_core_orbital(j,core_orb) ** 2 / rgrid(nt)%r(j)
      end do
      call simpson( size(work), work, rgrid(nt)%rab(:), integral )
      occupation = 2 * ( 2 * paw_recon(nt)%gipaw_core_orbital_l(core_orb) + 1 )
      nmr_shift_core(nt) = nmr_shift_core(nt) + occupation * integral
    enddo
    deallocate ( work )

    nmr_shift_core(nt) = nmr_shift_core(nt) * 17.75045395 * 1e-6
  enddo


  ! print integrals
  if (iverbosity > 10) then
    write(stdout,'(5X,''GIPAW integrals: -------------------------------------------'')')
    write(stdout,'(5X,''Atom  i/j   nmr_para   nmr_dia   epr_rmc  epr_para   epr_dia'')')
    do nt = 1, ntyp
      do il1 = 1, paw_recon(nt)%paw_nbeta
        l1 = paw_recon(nt)%psphi(il1)%label%l
        do il2 = 1, paw_recon(nt)%paw_nbeta
          l2 = paw_recon(nt)%psphi(il2)%label%l

          if (l1 /= l2) cycle
          if (il1 < il2) cycle

          write(stdout,1000) atm(nt), il1, il2, &
             radial_integral_paramagnetic(il1,il2,nt), &
             radial_integral_diamagnetic(il1,il2,nt), &
             radial_integral_rmc(il1,il2,nt), &
             radial_integral_paramagnetic_so(il1,il2,nt), &
             radial_integral_diamagnetic_so(il1,il2,nt)

        enddo
      enddo
    enddo
    write(stdout,'(5X,''------------------------------------------------------------'')')
    write(stdout,*)
  endif
1000 format(5X,A5,1X,I1,1X,I1,1X,5E10.2)

END SUBROUTINE gipaw_setup_integrals


!-----------------------------------------------------------------------
SUBROUTINE gipaw_setup_l
  !-----------------------------------------------------------------------
  !
  ! ... Setup the L operator using the properties of the cubic harmonics.
  ! ... Written by Ari P. Seitsonen and Uwe Gerstman
  !
  USE gipaw_module
  USE kinds,         ONLY : dp
  USE parameters,    ONLY : lmaxx
  USE io_global,     ONLY : stdout, ionode

  implicit none
  integer :: lm, l, m, lm1, lm2, m1, m2, abs_m1, abs_m2
  integer :: sign_m1, sign_m2
  real(dp) :: alpha_lm, beta_lm
  integer, allocatable :: lm2l(:),lm2m (:)
#ifdef DEBUG_CUBIC_HARMONIC
  real(dp) :: mysum1(3,lmaxx)
  real(dp) :: mysum2(3,1:lmaxx)
#endif

  ! L_x, L_y and L_z
  allocate ( lx(lmaxx**2,lmaxx**2) )
  allocate ( ly(lmaxx**2,lmaxx**2) )
  allocate ( lz(lmaxx**2,lmaxx**2) )
  allocate ( lm2l(lmaxx**2), lm2m(lmaxx**2) )
    
  lm = 0
  do l = 0, lmaxx - 1
    do m = 0, l
      lm = lm + 1
      lm2l ( lm ) = l
      lm2m ( lm ) = m
      if ( m /= 0 ) then
        lm = lm + 1
        lm2l ( lm ) = l
        lm2m ( lm ) = - m
      end if
    end do
  end do
    
  lx = 0.d0
  ly = 0.d0
  lz = 0.d0
  do lm2 = 1, lmaxx**2
    do lm1 = 1, lmaxx**2
      if ( lm2l ( lm1 ) /= lm2l ( lm2 ) ) cycle
          
      l = lm2l ( lm1 )
          
      m1 = lm2m ( lm1 )
      m2 = lm2m ( lm2 )
          
      ! L_x, L_y
      if ( m2 == 0 ) then
        if ( m1 == -1 ) then
          lx ( lm1, lm2 ) = - sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
        else if ( m1 == +1 ) then
          ly ( lm1, lm2 ) = + sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
        end if
      else if ( m1 == 0 ) then
        if ( m2 == -1 ) then
          lx ( lm1, lm2 ) = + sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
        else if ( m2 == +1 ) then
          ly ( lm1, lm2 ) = - sqrt(real(l*(l+1),dp)) / sqrt(2.0_dp)
        end if
      else
        abs_m1 = abs ( m1 )
        abs_m2 = abs ( m2 )
        sign_m1 = sign ( 1, m1 )
        sign_m2 = sign ( 1, m2 )
        alpha_lm = sqrt(real(l*(l+1)-abs_m2*(abs_m2+1),dp))
        beta_lm  = sqrt(real(l*(l+1)-abs_m2*(abs_m2-1),dp))
        if ( abs_m1 == abs_m2 + 1 ) then
          lx ( lm1, lm2 ) =-( sign_m2 - sign_m1 ) * alpha_lm / 4.0_dp
          ly ( lm1, lm2 ) = ( sign_m2 + sign_m1 ) * alpha_lm / 4.0_dp / sign_m2
        else if ( abs_m1 == abs_m2 - 1 ) then
          lx ( lm1, lm2 ) =-( sign_m2 - sign_m1 ) * beta_lm / 4.0_dp
          ly ( lm1, lm2 ) =-( sign_m2 + sign_m1 ) * beta_lm / 4.0_dp / sign_m2
        end if
      end if
          
      ! L_z
      if ( m1 == - m2 ) then
        lz ( lm1, lm2 ) = - m2
      end if
          
    end do
  end do
    
#ifdef DEBUG_CUBIC_HARMONICS
  write(stdout,'(A)') "lx:"
  write(stdout,'(9F8.5)') lx
       
  write(stdout,'(A)') "ly:"
  write(stdout,'(9F8.5)') ly
       
  write(stdout,'(A)') "lz:"
  write(stdout,'(9F8.5)') lz
       
  ! checks
  mysum1 = 0
  mysum2 = 0
  do lm2 = 1, lmaxx**2
    do lm1 = 1, lmaxx**2
      if ( lm2l ( lm1 ) /= lm2l ( lm2 ) ) cycle
      l = lm2l ( lm2 )
      mysum1(1,l+1) = mysum1(1,l+1) + lx(lm1,lm2)
      mysum2(1,l+1) = mysum2(1,l+1) + lx(lm1,lm2)**2
      mysum1(2,l+1) = mysum1(2,l+1) + ly(lm1,lm2)
      mysum2(2,l+1) = mysum2(2,l+1) + ly(lm1,lm2)**2
      mysum1(3,l+1) = mysum1(3,l+1) + lz(lm1,lm2)
      mysum2(3,l+1) = mysum2(3,l+1) + lz(lm1,lm2)**2
    end do
  end do
  write(stdout,'(A,9F8.4)') "Debug, sum1: x = ", mysum1(1,:)
  write(stdout,'(A,9F8.4)') "Debug, sum1: y = ", mysum1(2,:)
  write(stdout,'(A,9F8.4)') "Debug, sum1: z = ", mysum1(3,:)
  write(stdout,'(A,9F8.4)') "Debug, sum2: x = ", mysum2(1,:)
  write(stdout,'(A,9F8.4)') "Debug, sum2: y = ", mysum2(2,:)
  write(stdout,'(A,9F8.4)') "Debug, sum2: z = ", mysum2(3,:)
#endif
    
 deallocate ( lm2l, lm2m )

END SUBROUTINE gipaw_setup_l
    

