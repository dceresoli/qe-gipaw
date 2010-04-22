!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!-----------------------------------------------------------------------
SUBROUTINE hfi_fc_core_relax(fc_core)
  !-----------------------------------------------------------------------
  !
  ! ... Compute the core contribution to hyperfine Fermi contact
  !
  USE kinds,                 ONLY : dp
  USE constants,             ONLY : pi, tpi, fpi
  USE cell_base,             ONLY : tpiba, tpiba2
  USE parameters,            ONLY : ntypx
  USE ions_base,             ONLY : ntyp => nsp, atm, nat, tau, ityp, zv
  USE atom,                  ONLY : rgrid
  USE radial_grids,          ONLY : ndmx
  USE paw_gipaw,             ONLY : paw_recon
  USE scf,                   ONLY : rho
  USE gvect,                 ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ecutwfc
  USE gvect,                 ONLY : g, nl, gstart, ngm
  USE lsda_mod,              ONLY : nspin, isk, current_spin
  USE buffers,               ONLY : get_buffer
  USE gipaw_module,          ONLY : iverbosity
  USE klist,                 ONLY : nks, xk, wk
  USE wvfct,                 ONLY : npwx, nbnd, npw, igk, g2kin, wg, current_k
  USE becmod,                ONLY : calbec
  USE wavefunctions_module,  ONLY : evc
  USE io_global,             ONLY : stdout
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE mp_global,             ONLY : intra_pool_comm
  USE mp,                    ONLY : mp_sum
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE funct,                 ONLY : xc, xc_spin

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: fc_core(nat)

  ! -- constants ---------------------------------------------------------
  real(dp), parameter :: r_max = 10.d0      ! max atom radius
  integer, parameter :: n_max = 10          ! max number of s orbitals

  !-- local variables ----------------------------------------------------
  real(dp) :: eigenvalue(n_max,ntypx)
  real(dp) :: ae_orb(ndmx,n_max,ntypx)
  real(dp), allocatable :: vpot(:)
  integer :: nt, nn, zz, nstop
  real(dp) :: dx
  integer, external :: atomic_number

  integer :: s_maj, s_min, na, ispin, j
  complex(dp), allocatable :: aux(:)
  real(dp), allocatable :: sph_rho_bare(:,:,:), sph_rho_gipaw(:,:,:)

  integer :: il1, il2, nrc, l1, l2
  real(dp), allocatable :: rho_recon(:,:,:,:)

  integer :: ik, ibnd, ih, jh, ikb, jkb, m1, m2, lm1, lm2, ijkb0, nbs1, nbs2
  complex(dp) :: bec_product

  real(dp), allocatable :: delta_v(:,:), work(:)
  integer :: n1, n2, ncore, r_first, np
  real(dp) :: b(2), coeff, norm, contrib

  call start_clock('core_relax')
  fc_core = 0.d0

  !====================================================================
  ! recalculate AE orbitals (s only)
  !====================================================================
  write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
  if (iverbosity > 1) write(stdout,'(5X,''core-relax: calculating AE orbitals'')')
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

          ! solve atomic
          call ascheq(nn, 0, eigenvalue(nn,nt), rgrid(nt)%mesh, rgrid(nt), &
                      vpot, 2.d0*zz, 1d-12, ae_orb(1,nn,nt), nstop)
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

  !====================================================================
  ! project the density around each atom
  !====================================================================
  allocate( sph_rho_bare(ndmx,nat,nspin) )
  sph_rho_bare = 0.d0
  allocate( aux(nrxx) )
  do na = 1, nat
      nt = ityp(na)
      if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle
      if (iverbosity > 1) write(stdout,'(5X,''core-relax: projecting around atom '',I3)'), na
      ! bare density to reciprocal space
      do ispin = 1, nspin
        aux(1:nrxx) = rho%of_r(1:nrxx,ispin)
        call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
        call project_density(sph_rho_bare(:,na,ispin))
      enddo
  enddo
  call mp_sum(sph_rho_bare, intra_pool_comm)
  deallocate(aux)

  !====================================================================
  ! prepare reconstruction terms
  !====================================================================
  allocate( rho_recon(ndmx,paw_nkb,paw_nkb,ntyp) )
  rho_recon = 0.d0
  do nt = 1, ntyp
    !!if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle
    do il1 = 1, paw_recon(nt)%paw_nbeta
      nrc = paw_recon(nt)%psphi(il1)%label%nrc
      l1 = paw_recon(nt)%psphi(il1)%label%l
      if (l1 /= 0) cycle

      do il2 = 1, paw_recon(nt)%paw_nbeta
        l2 = paw_recon(nt)%psphi(il2)%label%l
        if (l2 /= 0 .or. l2 /= l1) cycle

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
  ! do GIPAW reconstruction
  !====================================================================
  allocate( sph_rho_gipaw(ndmx,nat,nspin) )
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
           do na = 1, nat
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
                       if (l2 /= 0 .or. l1 /= l2) cycle

                       bec_product = paw_becp(jkb,ibnd) &
                            * conjg( paw_becp(ikb,ibnd) )

                       sph_rho_gipaw(1:nrc,na,current_spin) = &
                            sph_rho_gipaw(1:nrc,na,current_spin) + &
                            rho_recon(1:nrc,nbs1,nbs2,nt) * &
                            bec_product * wg(ibnd,ik)
                    end do  ! jh
                 end do  ! ih
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              end if
           end do  ! na
        end do  ! nt
     end do  ! ibnd
  end do  ! ik

  !====================================================================
  ! compute perturbing potential
  !====================================================================
  allocate( delta_v(ndmx,nat) )
  delta_v = 0.d0
  do na = 1, nat
    nt = ityp(na)
    if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle
    do j = 1, rgrid(nt)%mesh
      if (rgrid(nt)%r(j) > r_max) cycle
      !! bare + gipaw
      b(1:2) = sph_rho_bare(j,na,1:2) + sph_rho_gipaw(j,na,1:2)
      delta_v(j,na) = -(1.d0/pi)*(b(s_maj)-b(s_min))/(b(1)+b(2))**(2.d0/3.d0)
      !! bare only
      !!b(1:2) = sph_rho_bare(j,na,1:2)
      !!delta_v(j,na) = -(2.d0/pi)*(b(s_maj)-b(s_min))/(b(1)+b(2))**(2.d0/3.d0)

      !! compute the true XC potential, but this overestimated core relaxation
      !!arho = abs(b(1)+b(2))
      !!zeta = (b(s_maj)-b(s_min))/arho
      !!call xc_spin(arho, zeta, ex, ec, vx(s_maj), vx(s_min), vc(s_maj), vc(s_min))
      !!delta_v(j,na) = vx(s_maj)+vc(s_maj)-vx(s_min)-vc(s_min)
      !!delta_v(j,na) = vx(s_maj)-vx(s_min)

      !! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
      !!write(70+na,'(5(F12.8))') rgrid(nt)%r(j), &
      !!   sph_rho_bare(j,na,1), sph_rho_bare(j,na,2), &
      !!   sph_rho_gipaw(j,na,1), sph_rho_gipaw(j,na,2)
      !!write(90+na,'(4(F12.8))') rgrid(nt)%r(j), delta_v(j,na), b(1), b(2)
      !! DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    enddo
  enddo

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

    do n1 = 1, ncore
      do n2 = n1+1, n_max
        if (eigenvalue(n2, nt) == 0.d0) cycle  ! unbound

        work = 0.d0
        do j = 1, rgrid(nt)%mesh
          work(j) = ae_orb(j,n1,nt) * delta_v(j,na) * ae_orb(j,n2,nt)
        end do
        call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab, coeff)
        contrib = 4.d0 * ae_orb(r_first,n1,nt) * ae_orb(r_first,n2,nt) * &
                  coeff / (eigenvalue(n1,nt) - eigenvalue(n2,nt)) / &
                  rgrid(nt)%r(r_first)**2   !/ fpi
        fc_core(na) = fc_core(na) + contrib
        if (iverbosity > 0) &
            write(stdout,'(5X,A,I3,2X,I1,''S -> '',I1,''S :'',F12.6)') &
                 atm(ityp(na)), na, n1, n2, contrib
      enddo
    enddo

  enddo  ! na

  write(stdout,*)

  deallocate(sph_rho_bare)
  deallocate(rho_recon)
  deallocate(sph_rho_gipaw)
  deallocate(delta_v)
  deallocate(work)

  call stop_clock('core_relax')

CONTAINS

  SUBROUTINE project_density(sph)
  implicit none
  real(dp), intent(inout) :: sph(:)
  complex(dp) :: rho0g
  real(dp) :: gg, gr, arg
  integer :: ig, jmax

  jmax = rgrid(nt)%mesh
  do j = 1, rgrid(nt)%mesh
      if (rgrid(nt)%r(j) < r_max) jmax = j
  enddo

  ! loop over g-vectors
  do ig = 1, ngm
    gg = sqrt(sum(g(:,ig)**2.d0))
    if (gg < 1.d-10) then  ! g = 0
      sph(:) = sph(:) + real(aux(nl(ig)))
    else
      arg = (tau(1,na)*g(1,ig)+tau(2,na)*g(2,ig)+tau(3,na)*g(3,ig))*tpi
      rho0g = aux(nl(ig)) * cmplx(cos(arg),sin(arg))

      do j = 1, jmax
        if (rgrid(nt)%r(j) < 1d-6) then
          sph(j) = sph(j) + real(rho0g)
        else
          gr = tpiba * gg * rgrid(nt)%r(j) 
          sph(j) = sph(j) + real(rho0g) * sin(gr)/gr
        endif
      enddo

    endif
  enddo
  END SUBROUTINE project_density

END SUBROUTINE hfi_fc_core_relax








#if 0
    ! loop over core orbitals
    ncore = paw_recon(nt)%gipaw_ncore_orbital
    do icore = 1, paw_recon(nt)%gipaw_ncore_orbital
      if (paw_recon(nt)%gipaw_core_orbital_l(icore) /= 0) cycle
      
      ! loop over core orbitals
      do icore2 = icore+1, paw_recon(nt)%gipaw_ncore_orbital
        if (paw_recon(nt)%gipaw_core_orbital_l(icore2) /= 0) cycle

        do j = 1, msh(nt)
          work(j) = paw_recon(nt)%gipaw_core_orbital(j,icore2) * &
                    delta_v(j,na) * &
                    paw_recon(nt)%gipaw_core_orbital(j,icore)
        end do
        call simpson(msh(nt), work, rab(:,nt), coeff)
        fc_core(na) = fc_core(na) + 4.d0 * &
          paw_recon(nt)%gipaw_core_orbital(r_first,icore) * &
          paw_recon(nt)%gipaw_core_orbital(r_first,icore2) * &
          coeff / (eigenvalues(icore,nt) - eigenvalues(icore2,nt)) / &
          r(r_first,nt)**2 / fpi
        write(stdout,'(5X,''CORE-CORE:'',2I3,10x,''fc_core up to now='',F12.6)')&
           icore, icore2, fc_core(na)
      enddo  ! ivale
    enddo  ! icore

    ! loop over core orbitals
    do icore = 1, paw_recon(nt)%gipaw_ncore_orbital
      if (paw_recon(nt)%gipaw_core_orbital_l(icore) /= 0) cycle
      
      ! loop over valence orbitals
      do ivale = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(ivale)%label%nrc
        l1 = paw_recon(nt)%psphi(ivale)%label%l
        if (l1 /= 0) cycle

        do j = 1, msh(nt)
          work(j) = paw_recon(nt)%aephi(ivale)%psi(j) * &
                    delta_v(j,na) * &
                    paw_recon(nt)%gipaw_core_orbital(j,icore)
        end do
        !call simpson(nrc, work, rab(:,nt), coeff)
        call simpson(msh(nt), work, rab(:,nt), coeff)
        !!!print*, eigenvalues(icore,nt), eigenvalues(ncore+ivale,nt)
        fc_core(na) = fc_core(na) + 4.d0 * &
          paw_recon(nt)%gipaw_core_orbital(r_first,icore) * &
          paw_recon(nt)%aephi(ivale)%psi(r_first) * &
          coeff / (eigenvalues(icore,nt) - eigenvalues(ncore+ivale,nt)) / &
          r(r_first,nt)**2 / fpi
        write(stdout,'(5X,''CORE-VALE:'',2I3,10x,''fc_core up to now='',F12.6)')&
           icore, ivale, fc_core(na)
      enddo  ! ivale
    enddo  ! icore
#endif

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

