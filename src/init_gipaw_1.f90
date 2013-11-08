!
! Copyright (C) 2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#undef DEBUG_THIS_ROUTINE

!----------------------------------------------------------------------
SUBROUTINE init_gipaw_1
  !----------------------------------------------------------------------
  !
  ! This routine initialize the variables of the paw projector
  ! and create the projectors in radial part (paw_beta) 
  !
  USE kinds ,      ONLY : dp
  USE parameters , ONLY : lqmax , lmaxx
  USE gipaw_module,ONLY : nbrx, pawproj
  USE cell_base ,  ONLY : omega
  USE ions_base,   ONLY : nat, ntyp => nsp, ityp, atm
  USE constants,   ONLY : fpi
  USE us,          ONLY : dq, nqx, tab, tab_d2y, qrad, spline_ps
  USE paw_gipaw,   ONLY : paw_recon, paw_nkb, paw_lmaxkb
  USE splinelib
  USE uspp,        ONLY : ap, aainit
  USE atom,        ONLY : rgrid, msh
  USE io_global,   ONLY : stdout
#ifdef __BANDS
  USE mp_bands,    ONLY : intra_bgrp_comm, inter_bgrp_comm
#endif
  USE mp_pools,    ONLY : intra_pool_comm
  USE mp,          ONLY : mp_sum
  USE uspp_param,  ONLY : upf  
  !--------------------------------------------------------------------
  implicit none
  integer :: nt, ih, jh, nb, l, m, ir, iq, startq
  integer :: lastq, na, j, n1, n2, ndm, nrs, nrc, lmaxkb
  real(dp), allocatable :: aux (:),aux1(:), besr(:)
  real(dp), allocatable :: s(:,:), sinv(:,:)
  real(DP), allocatable :: xdata(:)
  real(dp) :: prefr, pref, qi, norm
  real(dp) ::  vqint, rc, rs, pow, d1
  
  ! Initialization
  call start_clock ('init_gipaw_1')
  prefr = fpi / omega

  ndm = maxval(msh(1:ntyp))
  allocate (aux(ndm))
  allocate (aux1(ndm))
  allocate (besr(ndm))

  paw_lmaxkb = 0
  do nt = 1, ntyp
     lmaxkb = 0
     paw_recon(nt)%paw_nh = 0
     do nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        paw_recon(nt)%paw_nh = paw_recon(nt)%paw_nh + 2 * l + 1
        lmaxkb = max ( lmaxkb, l )
     enddo
     paw_lmaxkb = max ( paw_lmaxkb, lmaxkb )
     
     allocate (paw_recon(nt)%paw_nhtol(paw_recon(nt)%paw_nh))
     allocate (paw_recon(nt)%paw_nhtom(paw_recon(nt)%paw_nh))
     allocate (paw_recon(nt)%paw_indv(paw_recon(nt)%paw_nh))
     allocate (paw_recon(nt)%paw_tab(nqx,nbrx))
     allocate (paw_recon(nt)%paw_nl(0:lmaxkb))
     allocate (paw_recon(nt)%paw_iltonh(0:lmaxkb,paw_recon(nt)%paw_nh))
  enddo
  
  ! Calculate the number of beta functions of the solid
  paw_nkb = 0
  do na = 1, nat
     nt = ityp(na)
     paw_nkb = paw_nkb + paw_recon(nt)%paw_nh
  end do
  
  ! For each pseudopotential we initialize the indices nhtol, nhtom, indv
  write(stdout,'(5X,''GIPAW projectors -----------------------------------------------'')')
  do nt = 1, ntyp
     if (nt > 1) write(stdout,*)
     paw_recon(nt)%paw_nl = 0
     paw_recon(nt)%paw_iltonh = 0
     ih = 1
     do nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        paw_recon(nt)%paw_nl(l) = paw_recon(nt)%paw_nl(l) + 1
        paw_recon(nt)%paw_iltonh(l,paw_recon(nt)%paw_nl(l)) = nb
        do m = 1, 2 * l + 1
           paw_recon(nt)%paw_nhtol(ih) = l
           paw_recon(nt)%paw_nhtom(ih) = m
           paw_recon(nt)%paw_indv(ih) = nb
           ih = ih + 1
        end do
     end do
     
     ! Rescale the wavefunctions so that int_0^rc f|psi|^2=1
     pow = 1.0_dp
     do j = 1, paw_recon(nt)%paw_nbeta
        rc = paw_recon(nt)%psphi(j)%label%rc
        rs = 2.0_dp / 3.0_dp * rc
        nrc = COUNT ( rgrid(nt)%r(1:msh(nt)) <= rc )
        nrs = COUNT ( rgrid(nt)%r(1:msh(nt)) <= rs )
        write(stdout,'(5X,''atom='',A5,''  l='',I1,''  rc='',F10.4,''  rs='',F10.4)') &
              atm(nt), paw_recon(nt)%aephi(j)%label%l, rc, rs

        if (nrc < 1 .or. nrc > msh(nt)) call errore("init_gipaw_1", "impossible value for nrc", 1)
        if (nrs < 1 .or. nrs > msh(nt)) call errore("init_gipaw_1", "impossible value for nrs", 1)
        paw_recon(nt)%psphi(j)%label%nrc = nrc
        paw_recon(nt)%aephi(j)%label%nrc = nrc
        paw_recon(nt)%psphi(j)%label%nrs = nrs
        paw_recon(nt)%aephi(j)%label%nrs = nrs
        call step_f(aux, paw_recon(nt)%psphi(j)%psi**2, rgrid(nt)%r(:), &
             nrs, nrc, pow, msh(nt))
        call simpson(msh(nt), aux, rgrid(nt)%rab, norm)
        
        paw_recon(nt)%psphi(j)%psi = paw_recon(nt)%psphi(j)%psi / sqrt(norm)
        paw_recon(nt)%aephi(j)%psi = paw_recon(nt)%aephi(j)%psi / sqrt(norm)
     end do  ! j
     
     ! calculate the overlap matrix
     aux = 0.d0
     do l = 0, ubound(paw_recon(nt)%paw_nl,1)
        if (paw_recon(nt)%paw_nl(l) > 0) then

           allocate(s(paw_recon(nt)%paw_nl(l),paw_recon(nt)%paw_nl(l)))
           allocate(sinv(paw_recon(nt)%paw_nl(l),paw_recon(nt)%paw_nl(l)))
           do ih = 1, paw_recon(nt)%paw_nl(l)
              n1 = paw_recon(nt)%paw_iltonh(l,ih)
              do jh = 1, paw_recon(nt)%paw_nl(l)
                 n2 = paw_recon(nt)%paw_iltonh(l,jh)
                 
                 nrc = min(paw_recon(nt)%psphi(n1)%label%nrc, paw_recon(nt)%psphi(n2)%label%nrc)
                 nrs = min(paw_recon(nt)%psphi(n1)%label%nrs, paw_recon(nt)%psphi(n2)%label%nrs)
                 
                 call step_f(aux, paw_recon(nt)%psphi(n1)%psi(1:msh(nt)) &
                      * paw_recon(nt)%psphi(n2)%psi(1:msh(nt)), &
                      rgrid(nt)%r(:), nrs, nrc, pow, msh(nt))
                 
                 call simpson(msh(nt), aux, rgrid(nt)%rab, s(ih,jh))
                 
                 ! check if linearly dependent
                 if (ih < jh) then
                    if ( abs(abs(s(ih,jh)) - 1d0) < 1d-5 ) then
                       write(stdout,1000) l, ih, jh, s(ih,jh)
                       call flush_unit(stdout)
                       call errore("init_gipaw_1", "two projectors are linearly dependent", 1)
                    elseif ( abs(abs(s(ih,jh)) - 1d0) < 1d-2) then
                       write(stdout,1001) l, ih, jh, s(ih,jh)
                       call flush_unit(stdout)
                    endif
                 endif

              enddo  ! jh
           enddo  ! ih
1000 format(5X,'projs linearly dependent: l=',I1,'  n1,n2=',I2,',',I2,'  s=',F12.8)
1001 format(5X,'projs nearly linearly dependent: l=',I1,'  n1,n2=',I2,',',I2,'  s=',F12.8)
           
#ifdef DEBUG_THIS_ROUTINE
           !<apsi>
           if (iverbosity > 20) then
              do ih = 1, paw_recon(nt)%paw_nl(l)
                 do jh = ih, paw_recon(nt)%paw_nl(l)
                    write( stdout, '( A, I3, 3I2, F12.7 )' ) &
                         "PROJ: ", nt, l, ih, jh, s(ih,jh)
                 end do
              end do
           end if
           !</apsi>
#endif
           
           call invmat(paw_recon(nt)%paw_nl(l), s, sinv, norm)
           
           do ih = 1, paw_recon(nt)%paw_nl(l)
              n1 = paw_recon(nt)%paw_iltonh(l,ih)
              paw_recon(nt)%paw_betar(1:msh(nt),n1) = 0.d0
              
              do jh = 1, paw_recon(nt)%paw_nl(l)
                 n2 = paw_recon(nt)%paw_iltonh(l,jh)
                 paw_recon(nt)%paw_betar(1:msh(nt),n1) &
                      = paw_recon(nt)%paw_betar(1:msh(nt),n1) &
                      + sinv(ih,jh) * paw_recon(nt)%psphi(n2)%psi(1:msh(nt))
              enddo  ! jh
              
              nrc = paw_recon(nt)%psphi(n1)%label%nrc
              nrs = paw_recon(nt)%psphi(n1)%label%nrs
              
              call step_f(aux, paw_recon(nt)%paw_betar(1:msh(nt),n1), &
                   rgrid(nt)%r(:), nrs,nrc,pow,msh(nt))

              paw_recon(nt)%paw_betar(1:ndm,n1) = aux(1:ndm)

              ! PAW part (by Emine)
              if (pawproj(nt)) then
                 write(stdout,'(5X,''resetting projs and partial waves to PAW ones for atom '',A6)') atm(nt)
                 paw_recon(nt)%paw_betar(1:ndm,n1) = upf(nt)%beta(1:ndm,n1)
                 do nb = 1,paw_recon(nt)%paw_nbeta
                    paw_recon(nt)%aephi(nb)%psi(:) = upf(nt)%aewfc(:, nb)
                    paw_recon(nt)%psphi(nb)%psi(:) = upf(nt)%pswfc(:, nb)
                 enddo
              endif
           enddo  ! ih

           deallocate(sinv)
           deallocate(s)
        endif
     enddo  ! l

  enddo ! nt
  write(stdout,'(5X,''-----------------------------------------------------------------'')')
  write(stdout,*)

  ! compute Clebsch-Gordan coefficients
  call aainit(lmaxx+1)
  
  ! fill the interpolation table tab
  pref = fpi / sqrt ( omega )
  call divide (intra_pool_comm, nqx, startq, lastq)
  do nt = 1, ntyp
     paw_recon(nt)%paw_tab (:,:) = 0.d0
     do nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        do iq = startq, lastq
           qi = ( iq - 1 ) * dq
           call sph_bes ( msh(nt), rgrid(nt)%r, qi, l, besr )
           do ir = 1, msh(nt)
              aux(ir) = paw_recon(nt)%paw_betar(ir,nb) &
                   * besr(ir) * rgrid(nt)%r(ir)
           enddo
           call simpson ( msh(nt), aux, rgrid(nt)%rab, vqint )
           paw_recon(nt)%paw_tab(iq,nb) = vqint * pref
        enddo
     enddo
  
#ifdef __MPI
#  ifdef __BANDS
     call mp_sum ( paw_recon(nt)%paw_tab(:,:), intra_bgrp_comm )
     call mp_sum ( paw_recon(nt)%paw_tab(:,:), inter_bgrp_comm )
#  else
     call mp_sum ( paw_recon(nt)%paw_tab(:,:), intra_pool_comm )
#  endif
#endif

  enddo
  
  ! initialize spline interpolation
  if (spline_ps) then
     allocate (xdata(nqx))
     do iq = 1, nqx
        xdata(iq) = (iq - 1) * dq
     end do
     do nt = 1, ntyp
        allocate (paw_recon(nt)%paw_tab_d2y(nqx,paw_recon(nt)%paw_nbeta))
        paw_recon(nt)%paw_tab_d2y = 0.d0
        do nb = 1, paw_recon(nt)%paw_nbeta
           l = paw_recon(nt)%aephi(nb)%label%l
           d1 = ( paw_recon(nt)%paw_tab(2,nb) - paw_recon(nt)%paw_tab(1,nb) ) &
                / dq
           call spline(xdata, paw_recon(nt)%paw_tab(:,nb), 0.d0, d1, &
                paw_recon(nt)%paw_tab_d2y(:,nb))
        end do
     end do
     deallocate (xdata)
  end if
  
  deallocate (besr)
  deallocate (aux1)
  deallocate (aux)
  
  call stop_clock ('init_gipaw_1')
  return
  
END SUBROUTINE init_gipaw_1




!----------------------------------------------------------------------
SUBROUTINE step_f(f2,f,r,nrs,nrc,pow,mesh)
  !----------------------------------------------------------------------
  !
  ! This routine apply a function which go smoothly to zero from rs to rc
  ! 
  USE kinds, ONLY : dp
  implicit none
  integer, intent(in) :: mesh
  real(dp), intent(out):: f2(mesh)
  real(dp), intent(in) :: f(mesh), r(mesh), pow
  real(dp) :: rcp, rsp
  integer :: nrs, nrc 
  integer :: i
 
  rcp = r(nrc)
  rsp = r(nrs)

  do i = 1, mesh
     if (r(i) <= rsp) then
        f2(i) = f(i)
     elseif (r(i) <= rcp) then
        f2(i) = f(i) * (1.d0-3.d0*((r(i)-rsp)/(rcp-rsp))**2 + &
                        2.d0*((r(i)-rsp)/(rcp-rsp))**3)**pow
     else
        f2(i) = 0.d0
     endif
  enddo

END SUBROUTINE step_f

