!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!-----------------------------------------------------------------------
SUBROUTINE hfi_fc_core_relax(method, fc_core)
  !-----------------------------------------------------------------------
  !
  ! ... Compute the core contribution to hyperfine Fermi contact
  ! ... according to: Bahramy, Sluiter, Kawazoe, PRB 76, 035124 (2007)
  !
  USE kinds,                 ONLY : dp
  USE constants,             ONLY : pi, tpi, fpi, e2
  USE cell_base,             ONLY : tpiba, tpiba2
  USE parameters,            ONLY : ntypx
  USE ions_base,             ONLY : ntyp => nsp, atm, nat, tau, ityp, zv
  USE atom,                  ONLY : rgrid
  USE radial_grids,          ONLY : ndmx
  USE paw_gipaw,             ONLY : paw_recon
  USE scf,                   ONLY : rho
  USE gvect,                 ONLY : g, nl, gstart, ngm
  USE fft_base,              ONLY : dfftp
  USE fft_interfaces,        ONLY : fwfft
  USE lsda_mod,              ONLY : nspin, isk, current_spin
  USE buffers,               ONLY : get_buffer
  USE gipaw_module,          ONLY : iverbosity
  USE klist,                 ONLY : nks, xk, wk
  USE wvfct,                 ONLY : npwx, nbnd, npw, igk, g2kin, wg, current_k, ecutwfc
  USE becmod,                ONLY : calbec
  USE wavefunctions_module,  ONLY : evc
  USE io_global,             ONLY : stdout
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE mp_pools,              ONLY : intra_pool_comm, inter_pool_comm
  USE mp,                    ONLY : mp_sum
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE funct,                 ONLY : xc, xc_spin
  USE gipaw_module,          ONLY : nbrx

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  integer, intent(in) :: method
  real(dp), intent(out) :: fc_core(nat)

  ! -- constants ---------------------------------------------------------
  !real(dp), parameter :: r_max = 2.5d0      ! max core radius
  real(dp), parameter :: r_max = 5.0d0      ! max core radius
  integer, parameter :: n_max = 10          ! max number of s orbitals

  !-- local variables ----------------------------------------------------
  real(dp) :: eigenvalue(n_max,ntypx)
  real(dp) :: ae_orb(ndmx,n_max,ntypx)
  real(dp), allocatable :: vpot(:)
  integer :: nt, nn, zz, nstop
  integer, external :: atomic_number

  integer :: s_maj, s_min, na, ispin, j
  complex(dp), allocatable :: aux(:), rho_g(:)
  real(dp), allocatable :: sph_rho_bare(:,:), sph_rho_gipaw(:,:)
  real(dp) :: sph_rho_core(ndmx)

  integer :: il1, il2, nrc, l1, l2
  real(dp), allocatable :: rho_recon(:,:,:,:)

  integer :: ik, ibnd, ih, jh, ikb, jkb, m1, m2, lm1, lm2, ijkb0, nbs1, nbs2
  complex(dp) :: bec_product

  real(dp), allocatable :: delta_v(:,:), work(:)
  integer :: n1, n2, ncore, r_first
  real(dp) :: b(2), coeff, norm, contrib
  integer :: mode, nin, mesh
  real(dp) :: arho, zeta, ex, ec, vx(2), vc(2)

  real(dp), allocatable :: rho_(:,:), rhoc(:), vgc(:,:), egc(:,:), tau_(:,:), vtau(:)

  if (method < 1 .or. method > 3) call errore('core-relax', 'unknown method', abs(method))

  call start_clock('core_relax')
  fc_core = 0.d0

  !====================================================================
  ! recalculate AE orbitals (s only)
  !====================================================================
  write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
  if (iverbosity > 1) write(stdout,'(5X,''core-relax: calculating AE orbitals'')')
  if (iverbosity > 1) write(stdout,'(5X,''core-relax: method = '',I1)') method
  eigenvalue(:,:) = 0.d0
  ae_orb(:,:,:) = 0.d0
  do nt = 1, ntyp
      if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle
      allocate( vpot(rgrid(nt)%mesh), work(rgrid(nt)%mesh) )
      do nn = 1, n_max
          ! setup grid
          zz = atomic_number(atm(nt))
          rgrid(nt)%dx = log( rgrid(nt)%r(2)/rgrid(nt)%r(1) )
          rgrid(nt)%r2 = rgrid(nt)%r**2
          rgrid(nt)%sqr = sqrt(rgrid(nt)%r)

          ! setup potential
          vpot(:) = paw_recon(nt)%gipaw_ae_vloc(:) / rgrid(nt)%r(:)

          ! solve radial schroedinger equation
          !!call ascheq(nn, 0, eigenvalue(nn,nt), rgrid(nt)%mesh, rgrid(nt), &
          !!            vpot, 2.d0*zz, 1d-12, ae_orb(1,nn,nt), nstop)
          mode = 2                ! non-relativistic
          if (zz >= 20) mode = 1  ! scalar-relativistic
          eigenvalue(nn,nt) = -(dble(zz)/dble(nn))**2.0
          nin = 0
          call lschps(mode, 2.d0*zz, 1d-12, rgrid(nt), nin, nn, 0, &
                      eigenvalue(nn,nt), vpot, ae_orb(1,nn,nt), nstop)
          if (nstop /= 0 .and. nstop /= 50) then
              eigenvalue(nn,nt) = 0.d0
              exit
          endif 

          if (iverbosity > 0) then
              do j = 1, rgrid(nt)%mesh
                  work(j) = ae_orb(j,nn,nt)**2.d0
                  !!write(90+nn,*) rgrid(nt)%r(j), ae_orb(j,nn,nt)
              enddo
              call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab, norm)
              write(stdout,'(5X,A,4X,I1,''S   eig='',F12.4,'' Ry   norm='',F12.4)') &
                    atm(nt), nn, eigenvalue(nn,nt), norm
          endif

      enddo
      deallocate( vpot, work )
      if (iverbosity > 1) write(stdout,*)
  enddo

  ! Select majority and minority spin components
  call select_spin(s_min, s_maj)

  ! prepare reconstruction terms
  allocate( rho_recon(ndmx,nbrx,nbrx,ntyp) )
  rho_recon = 0.d0
  do nt = 1, ntyp

    do il1 = 1, paw_recon(nt)%paw_nbeta
      nrc = paw_recon(nt)%psphi(il1)%label%nrc
      l1 = paw_recon(nt)%psphi(il1)%label%l
      if (l1 /= 0) cycle

      do il2 = 1, paw_recon(nt)%paw_nbeta
        l2 = paw_recon(nt)%psphi(il2)%label%l
        if (l2 /= 0) cycle

        do j = 1, nrc
          rho_recon(j,il1,il2,nt) = &
                   ( paw_recon(nt)%aephi(il1)%psi(j) &
                   * paw_recon(nt)%aephi(il2)%psi(j) &
                   - paw_recon(nt)%psphi(il1)%psi(j) &
                   * paw_recon(nt)%psphi(il2)%psi(j) ) &
                   / rgrid(nt)%r(j) ** 2 / fpi
        enddo

      enddo ! il2
    enddo ! il1
   enddo ! nt


  !====================================================================
  ! loop over atoms with core electrons
  !====================================================================
  allocate( sph_rho_bare(ndmx,nspin) )
  allocate( sph_rho_gipaw(ndmx,nspin) )
  allocate( aux(dfftp%nnr), rho_g(ngm) )
  allocate( delta_v(ndmx,nat) )
  delta_v = 0.d0

  if (iverbosity > 1) write(stdout,'(5X,''core-relax: max core radius (r_max) =  '',F12.4)') r_max
  do na = 1, nat
    nt = ityp(na)
    if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle
    if (iverbosity > 1) write(stdout,'(5X,''core-relax: projecting around atom '',I3)') na

    !====================================================================
    ! project the density around each atom
    !====================================================================
    sph_rho_bare = 0.d0
    do ispin = 1, nspin
      aux(1:dfftp%nnr) = rho%of_r(1:dfftp%nnr,ispin)
      call fwfft('Dense',aux,dfftp)
      rho_g(1:ngm) = aux(nl(1:ngm))
      call spherical_average(rgrid(nt)%mesh, rgrid(nt)%r, tau(1,na), r_max, rho_g, sph_rho_bare(1,ispin))
    enddo
    call mp_sum(sph_rho_bare, intra_pool_comm)
    ! the density coming out are such that: \int_0^\infty 4\pi r^2 \rho(r) dr = N_{el}

    !====================================================================
    ! do GIPAW reconstruction
    !====================================================================
    sph_rho_gipaw = 0.d0
    do ik = 1, nks
      current_k = ik
      current_spin = isk(ik)
     
      call gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      call get_buffer (evc, nwordwfc, iunwfc, ik)
      call init_gipaw_2 (npw, igk, xk(1,ik), paw_vkb)
      call calbec (npw, paw_vkb, evc, paw_becp)
     
      do ibnd = 1, nbnd
        ijkb0 = 0
        do nt = 1, ntyp
              if ( ityp(na) == nt ) then
                 do ih = 1, paw_recon(nt)%paw_nh
                    ikb = ijkb0 + ih
                    nbs1 = paw_recon(nt)%paw_indv(ih)
                    l1 = paw_recon(nt)%paw_nhtol(ih)
                    m1 = paw_recon(nt)%paw_nhtom(ih)
                    lm1 = m1 + l1**2
                    nrc = paw_recon(nt)%psphi(nbs1)%label%nrc
                    if (l1 /= 0) cycle
 
                    do jh = 1, paw_recon(nt)%paw_nh
                       jkb = ijkb0 + jh
                       nbs2 = paw_recon(nt)%paw_indv(jh)
                       l2 = paw_recon(nt)%paw_nhtol(jh)
                       m2 = paw_recon(nt)%paw_nhtom(jh)
                       lm2 = m2 + l2**2 
                       if (l2 /= 0) cycle

                       bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))

                       sph_rho_gipaw(1:nrc,current_spin) = &
                            sph_rho_gipaw(1:nrc,current_spin) + &
                            rho_recon(1:nrc,nbs1,nbs2,nt) * &
                            bec_product * wg(ibnd,ik)
                    end do  ! jh
                 end do  ! ih
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              end if
        end do  ! nt
      end do  ! ibnd
    end do  ! ik
    call mp_sum(sph_rho_gipaw, inter_pool_comm)

    !====================================================================
    ! compute perturbing potential
    !====================================================================
    nt = ityp(na)
    if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle

    mesh = rgrid(nt)%mesh

    sph_rho_core = 0.d0
    do n1 = 1, paw_recon(nt)%gipaw_ncore_orbital
      if (paw_recon(nt)%gipaw_core_orbital_l(n1) /= 0) cycle
      do j = 1, mesh
        sph_rho_core(j) = sph_rho_core(j) + 2.d0*paw_recon(nt)%gipaw_core_orbital(j,n1)**2.d0 / rgrid(nt)%r2(j) / fpi
      enddo
    enddo

    allocate(work(mesh))
    do j = 1, mesh
      work(j) = sph_rho_bare(j,1)+sph_rho_bare(j,2)+sph_rho_gipaw(j,1)+sph_rho_gipaw(j,2)+sph_rho_core(j)
      work(j) = work(j) * fpi * rgrid(nt)%r(j)**2.0

      if (rgrid(nt)%r(j) > r_max) cycle

      ! calculate density and magnetization
      b(1:2) = sph_rho_bare(j,1:2) + sph_rho_gipaw(j,1:2) + sph_rho_core(j)
      !b(1:2) = b(1:2) * rgrid(nt)%r(j)**2 * fpi
      arho = abs(b(1)+b(2))
      zeta = (b(s_maj)-b(s_min))/arho

      ! compute the perturbing potential, three possibilities
      select case (method)
         case (1) ! simple local exchange: Eq.(20) or PRB 76, 035124
         delta_v(j,na) = -(e2*2.d0/pi) * (b(s_maj)-b(s_min)) / arho**(2.d0/3.d0)

         case (2) ! Exchange only
         call xc_spin(arho, zeta, ex, ec, vx(s_maj), vx(s_min), vc(s_maj), vc(s_min))
         delta_v(j,na) = e2*(vx(s_maj) - vx(s_min))

         case (3) ! Full XC
         call xc_spin(arho, zeta, ex, ec, vx(s_maj), vx(s_min), vc(s_maj), vc(s_min))
         delta_v(j,na) = e2*(vx(s_maj)+vc(s_maj)-vx(s_min)-vc(s_min))
      end select

      !! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
      write(70+na,'(7(F12.6))') rgrid(nt)%r(j), sph_rho_bare(j,1), sph_rho_bare(j,2), &
                                sph_rho_gipaw(j,1), sph_rho_gipaw(j,2), sph_rho_core(j), delta_v(j,na)
      !! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    enddo

    ! gradient correction    
    if (method == 3) then
      allocate(rho_(mesh,2), rhoc(mesh), vgc(mesh,2), egc(mesh,2), tau_(mesh,2), vtau(mesh))
      rhoc = 0.d0
      tau_ = 0.d0
      vtau = 0.d0
      rho_(1:mesh,1) = (sph_rho_bare(1:mesh,1) + sph_rho_gipaw(1:mesh,1) + sph_rho_core(1:mesh)) * fpi * rgrid(nt)%r2(1:mesh)
      rho_(1:mesh,2) = (sph_rho_bare(1:mesh,2) + sph_rho_gipaw(1:mesh,2) + sph_rho_core(1:mesh)) * fpi * rgrid(nt)%r2(1:mesh)
      call vxcgc(mesh, mesh, nspin, rgrid(nt)%r, rgrid(nt)%r2, rho_, rhoc, vgc, egc, tau_, vtau, 1)
      delta_v(:,na) = delta_v(:,na) + vgc(:,1) - vgc(:,2)
      do j = 1, mesh
        write(80+na,'(4(F12.6))') rgrid(nt)%r(j), delta_v(j,na), vgc(j,1), vgc(j,2)
      enddo
      deallocate(rho_, rhoc, vgc, egc, tau_, vtau)
    endif

    call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab, norm)
    deallocate(work)
    if (iverbosity > 1) write(stdout,'(5X,''core-relax: integrated charge = '',F10.4)') norm

  !====================================================================
  ! end of the loop over atoms
  !====================================================================
  enddo
  deallocate(aux, rho_g, sph_rho_bare, sph_rho_gipaw)

  !====================================================================
  ! now, do the core relaxation via pertubation theory (PRB 76, 035124)
  !====================================================================
  allocate( work(ndmx) )
  fc_core(1:nat) = 0.d0

  do na = 1, nat
    nt = ityp(na)
    if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle

    do j = 1, rgrid(nt)%mesh
      if (rgrid(nt)%r(j) > 1d-5) then
        r_first = j
        exit
      endif
    enddo

    ! count number of s core orbitals
    ncore = 0
    do n1 = 1, paw_recon(nt)%gipaw_ncore_orbital
      if (paw_recon(nt)%gipaw_core_orbital_l(n1) /= 0) cycle
      ncore = ncore + 1
    enddo

    if (iverbosity > 1) then
      do n1 = 1, ncore
        work = 0.d0
        do j = 1, rgrid(nt)%mesh
          work(j) = ae_orb(j,n1,nt) * delta_v(j,na) * ae_orb(j,n1,nt)
        end do
        call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab, coeff)
        write(stdout,'(5X,A,I3,2X,I1''S splitting:'',F12.6)') atm(ityp(na)), na, n1, coeff
      enddo
    endif

    do n1 = 1, ncore
      do n2 = n1+1, n_max
        if (eigenvalue(n2, nt) == 0.d0) cycle  ! unbound

        work = 0.d0
        !!DEBUG: write(80+na,'(''#'',2I4)') n1, n2
        do j = 1, rgrid(nt)%mesh
          work(j) = ae_orb(j,n1,nt) * delta_v(j,na) * ae_orb(j,n2,nt)
          !!DEBUG: write(80+na,'(5(F12.8))') rgrid(nt)%r(j), ae_orb(j,n1,nt), ae_orb(j,n2,nt), &
          !!DEBUG:                           delta_v(j,na), work(j)
        end do
        !!DEBUG: write(80+na,*)

        call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab, coeff)
        contrib = 2.d0 * 4.d0 * ae_orb(r_first,n1,nt) * ae_orb(r_first,n2,nt) * &
                  coeff / (eigenvalue(n1,nt) - eigenvalue(n2,nt)) / &
                  rgrid(nt)%r(r_first)**2 / fpi

        fc_core(na) = fc_core(na) + contrib
        if (iverbosity > 0) &
            write(stdout,'(5X,A,I3,2X,I1,''S -> '',I1,''S :'',F12.6)') &
                 atm(ityp(na)), na, n1, n2, contrib
      enddo
    enddo

  enddo  ! na

  write(stdout,*)

  deallocate(rho_recon)
  deallocate(delta_v)
  deallocate(work)

  call stop_clock('core_relax')

END SUBROUTINE hfi_fc_core_relax





#if 0
    !--------------------------------------------------------------------
    ! this part recalculates the core orbitals
    !--------------------------------------------------------------------
    ! loop over core orbitals
    nn = 0
    do icore = 1, paw_recon(nt)%gipaw_ncore_orbital
      if (paw_recon(nt)%gipaw_core_orbital_l(icore) /= 0) cycle
      nn = nn + 1
      dx = log(r(2,nt)/r(1,nt))
      zz = get_z(atm(nt))

      ! averaged potential
      do j = 1, msh(nt)
        b(1) = 2.d0*(paw_recon(nt)%gipaw_core_orbital(j,icore)/r(j,nt))**2.0/fpi
        arho = abs(b(1))
        call xc(arho, ex, ec, vx(1), vc(1))
        vpot1(j) = vx(1) + vc(1)
        vpot1(j) = 0.d0
        !!write(92,*) r(j,nt), vpot1(j)
      enddo

      ! solve for majority spin
      vpot(:msh(nt)) = paw_recon(nt)%gipaw_ae_vloc(:msh(nt))/r(:msh(nt),nt) + &
                       0.5d0*delta_v(:msh(nt),na)
      call ascheq(nn, 0, ene, msh(nt), dx, r(:,nt), r(:,nt)**2, &
        sqrt(r(:,nt)), vpot, 2.d0*zz, 1d-6, y, nstop)
      print*, ene, zz, dx
      fc_core(na) = fc_core(na) + (y(r_first)/r(r_first,nt))**2 / fpi
      do j = 1, msh(nt)
        write(60+icore,*) r(j,nt), y(j)
      enddo

      ! solve for minority spin
      vpot(:msh(nt)) = paw_recon(nt)%gipaw_ae_vloc(:msh(nt))/r(:msh(nt),nt) - &
                       0.5d0*delta_v(:msh(nt),na)
      call ascheq(nn, 0, ene, msh(nt), dx, r(:,nt), r(:,nt)**2, &
        sqrt(r(:,nt)), vpot, zz, 1d-6, y, nstop)
      print*, ene
      fc_core(na) = fc_core(na) - (y(r_first)/r(r_first,nt))**2 / fpi
      do j = 1, msh(nt)
        write(65+icore,*) r(j,nt), y(j)
      enddo
    enddo
#endif

