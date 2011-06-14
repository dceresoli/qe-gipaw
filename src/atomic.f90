!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! slightly modified ascheq, lschps and vxcgc from ld1

!
!---------------------------------------------------------------
subroutine ascheq(nn,lam,e,mesh,grid,vpot,ze2,thresh0,y,nstop)
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation for
  !  bound states in a local potential.
  !  thresh determines the absolute accuracy for the eigenvalue
  !
  use kinds, only : DP
  use radial_grids, only: radial_grid_type, series
  implicit none
  type(radial_grid_type), intent(in) :: grid
  integer :: mesh,lam, ierr
  integer:: nn,nstop,maxter,iter,l1,i,ik=0,ncross,n, &
       nstart,ns,n2,nst2,ndcr
  real(DP) :: ze2,ddx12,eup,elw,ymx,rap,rstart,di,expn,  &
       c1,c2,c3,c4,fe,a0,a1,a2,sum0,f2,sum,sqlhf,f0,f1,dfe,de,eps,&
       yln,xp,sum1
  real(DP):: vpot(mesh), y(mesh)
  real(DP),allocatable:: c(:), el(:), f(:)
  real(DP):: b(0:3),e,thresh0, thresh
  data maxter/200/
  !
  !  set up constants and initialize
  !
  allocate(c(mesh),stat=ierr)
  allocate(f(mesh),stat=ierr)
  allocate(el(mesh),stat=ierr)

  thresh=thresh0
  !!if (e<-5.e+2) thresh=thresh0*10.0_DP
  iter=0
  ddx12=grid%dx*grid%dx/12.0_dp
  l1=lam+1
  sqlhf=(DBLE(lam)+0.5_dp)**2
  ndcr=nn-lam-1
  !
  !  set initial lower and upper bounds to the eigenvalue
  !
  eup=vpot(mesh)+sqlhf/grid%r2(mesh)
  elw=eup
  do i=1,mesh
     elw=min(elw,vpot(i)+sqlhf/grid%r2(i))
  enddo
  nstop=200
  if(eup.eq.elw) go to 900
  e = (eup + elw)/2.d0
  !!if(e.gt.eup) e=0.9_DP*eup+0.1_DP*elw
  !!if(e.lt.elw) e=0.9_DP*elw+0.1_DP*eup
  !
  !  series developement of the potential near the origin
  !
  do i=1,4
     y(i)=vpot(i)-ze2/grid%r(i)
  enddo
  call series(y,grid%r,grid%r2,b)
  !
300 continue
  iter=iter+1
  nstop=300
  if(iter.gt.maxter) go to 900
  !
  !  set up the f-function and determine the position of its last
  !  change of sign
  !  f < 0 (approximatively) means classically allowed   region
  !  f > 0         "           "        "      forbidden   "
  !
  f(1)=ddx12*(grid%r2(1)*(vpot(1)-e)+sqlhf)
  do i=2,mesh
     f(i)=ddx12*(grid%r2(i)*(vpot(i)-e)+sqlhf)
     if( f(i) .ne. sign(f(i),f(i-1)) ) ik=i
  enddo
  nstop=302
  
  if(ik.ge.mesh-2) go to 900
  do i=1,mesh
     f(i)=1.0_dp-f(i)
  enddo
  !
  y(:) = 0.0_dp
  !
  !  determination of the wave-function in the first two points by
  !  series developement
  !
  call start_scheq( lam, e, b, grid, ze2, y, c1, c2, c3, c4 )
  !
  !  start outward integration and count number of crossings
  !
  ncross=0
  ymx=0.0_dp
  do n=2,ik-1
     y(n+1)=((12.0_dp-10.0_dp*f(n))*y(n)-f(n-1)*y(n-1))/f(n+1)
     if ( y(n) .ne. sign(y(n),y(n+1)) ) ncross=ncross+1
     ymx=max(ymx,abs(y(n+1)))
  end do
  !
  !  matching radius has been reached going out. if ncross is not
  !  equal to ndcr, modify the trial eigenvalue.
  !
  !!WRITE(6,'(''iter='',I3,2X,''ncross='',I2,4X,''e,eup,elw='',3(F12.6,2X))') iter, ncross, e, eup, elw
  if(ndcr < ncross) then
     !
     !  too many crossings. e is an upper bound to the true eigen-
     !  value. increase abs(e)
     !
     eup=e
     rap=(DBLE(ncross+l1)/DBLE(nn))**2
     e=(e-vpot(mesh))*rap+vpot(mesh)
     if(e.lt.elw) e=0.9_dp*elw+0.1_dp*eup
     go to 300
  else if (ndcr > ncross) then
     !
     !  too few crossings. e is a lower bound to the true eigen-
     !  value. decrease abs(e)
     !
     elw=e
     rap=(DBLE(ncross+l1)/DBLE(nn))**2
     e=(e-vpot(mesh))*rap+vpot(mesh)
     if(e.gt.eup) e=0.9_dp*eup+0.1_dp*elw
     go to 300
  end if
  !
  !  prepare inward integration
  !  charlotte froese can j phys 41,1895(1963)
  !
  !            start at  min( rmax, 10*rmatch )
  !
  nstart=mesh
  ns=10
  rstart=ns*grid%r(ik)
  if(rstart.lt.grid%r(mesh)) then
     do  i=ik,mesh
        nstart=i
        if(grid%r(i).ge.rstart) go to 403
     enddo
403  nstart=nstart/2
     nstart=2*nstart+1
  end if
  !
  !  set up a, l, and c vectors
  !
  n=ik+1
  el(n)=10.0_dp*f(n)-12.0_dp
  c(n)=-f(ik)*y(ik)
  n2=ik+2
  do n=n2,nstart
     di=10.0_dp*f(n)-12.0_dp
     el(n)=di-f(n)*f(n-1)/el(n-1)
     c(n)=-c(n-1)*f(n-1)/el(n-1)
  enddo
  !
  !  start inward integration by the froese's tail procedure
  !
  expn=exp(-sqrt(12.0_dp*abs(1.0_dp-f(nstart-1))))
  y(nstart-1)=c(nstart-1)/(el(nstart-1)+f(nstart)*expn)
  y(nstart)=expn*y(nstart-1)
  do n=nstart-2,ik+1,-1
    y(n)=(c(n)-f(n+1)*y(n+1))/el(n)
 end do
  !
  !  if necessary, improve the trial eigenvalue by the cooley's
  !  procedure. jw cooley math of comp 15,363(1961)
  !
  fe=(12.0_dp-10.0_dp*f(ik))*y(ik)-f(ik-1)*y(ik-1)-f(ik+1)*y(ik+1)
  !
  !  calculate the normalization
  !
  if(ymx.ge.1.0e10_dp) then
     do  i=1,mesh
        y(i)=y(i)/ymx
     enddo
  end if
  a0=1.0_dp/DBLE(2*lam+3)
  a1=c1/DBLE(lam+2)
  a2=(c1*c1+c2+c2)/DBLE(2*lam+5)
  sum0=(a0+grid%r(1)*(a1+grid%r(1)*a2))*grid%r(1)**(2*lam+3)
  nst2=nstart-2
  f2=grid%r2(1  )*y(1  )*y(1  )
  sum=grid%r(1)*f2/DBLE(2*l1+1)
  do n=1,nst2,2
     f0=f2
     f1=grid%r2(n+1)*y(n+1)*y(n+1)
     f2=grid%r2(n+2)*y(n+2)*y(n+2)
     sum=sum+f0+f2+4.0_DP*f1
  enddo
  sum=sum0+grid%dx*sum/3.0_dp
  dfe=-y(ik)*f(ik)/grid%dx/sum
  de=-fe*dfe
  eps=abs(de/e)
  !!WRITE(6,'(''iter='',I2,2X,''n,l='',2(I1,1X),4X,''e,de='',F12.6,2X,E12.4)') iter, nn, lam, e, de
  if(abs(de).lt.thresh) go to 600
  if(eps.gt.0.25_dp) de=0.25_dp*de/eps
  if(de.gt.0.0_dp) elw=e
  if(de.lt.0.0_dp) eup=e
  e=e+de
  if(e.gt.eup) e=0.9_dp*eup+0.1_dp*elw
  if(e.lt.elw) e=0.9_dp*elw+0.1_dp*eup
  if(iter.lt.maxter) go to 300
  nstop=50
600 continue
  !
  !  normalize the eigenfunction and exit
  !
  do n=nstart,mesh-1
     y(n+1)=0.0_dp
     if(y(n).eq.0.0_dp) go to 601
     yln=log(abs(y(n)))
     xp=-sqrt(12.0_dp*abs(1.0_dp-f(n)))
     expn=yln+xp
     if(expn.lt.-80.0_dp) go to 601
     y(n+1)=sign(exp(expn),y(n))
601 continue
  enddo
  sum1=0.0_dp
  do n=nstart,mesh-2,2
     f0=f2
     f1=grid%r2(n+1)*y(n+1)*y(n+1)
     f2=grid%r2(n+2)*y(n+2)*y(n+2)
     sum1=sum1+f0+f2+4.0_dp*f1
  enddo
  sum=sum+grid%dx*sum1/3.0_dp
  sum=sqrt(sum)
  do n=1,mesh
     y(n)=grid%sqr(n)*y(n)/sum
  enddo
  if(nstop.lt.100) go to 900
  nstop=0
  deallocate(el)
  deallocate(f )
  deallocate(c )
  return
  !
  !  error exit
  !
  ! 900  write(6,9000) nstop,nn,lam,elw,eup
  ! 9000 format(5x,'error in ascheq: nstop =',i4,'. n l =',2i3,/ &
  !     & 5x,'elw =',f15.10,' eup =',f15.10)
900 continue
  deallocate(el)
  deallocate(f )
  deallocate(c )
  return

end subroutine ascheq


!---------------------------------------------------------------
subroutine start_scheq(lam,e,b,grid,ze2,y, c1,c2,c3,c4)
  !---------------------------------------------------------------
  !
  !  determines the wave-function in the first two points by
  !  series developement. It receives as input:
  !  lam the angular momentum 
  !  e   the energy in Ry 
  !  b(0:3) the coefficients of a polynomial that interpolates the 
  !         potential in the first three points
  !  grid  the mesh
  !  ze2   the zed of the mesh
  !  in output y(1:2) contains the solution in the first two points
  !
USE kinds, ONLY : DP
USE radial_grids, ONLY : radial_grid_type
IMPLICIT NONE
TYPE(radial_grid_type), INTENT(IN) :: grid
INTEGER, INTENT(IN) :: lam
REAL(DP), INTENT(IN) :: b(0:3), e
REAL(DP), INTENT(OUT) :: c1, c2, c3, c4
REAL(DP) :: ze2, xl1, x4l6, x6l12, x8l20, b0e, rr1, rr2
REAL(DP) :: y(1:2)
INTEGER :: l1
!
!  set up constants and initialize
!
l1=lam+1
xl1=lam+1.0_DP
x4l6=4.0_dp*lam+6.0_dp
x6l12=6.0_dp*lam+12.0_dp
x8l20=8.0_dp*lam+20.0_dp
!
!
b0e=b(0)-e
c1=0.5_dp*ze2/xl1
c2=(c1*ze2+b0e)/x4l6
c3=(c2*ze2+c1*b0e+b(1))/x6l12
c4=(c3*ze2+c2*b0e+c1*b(1)+b(2))/x8l20
rr1=(1.0_dp+grid%r(1)*(c1+grid%r(1)*(c2+grid%r(1)*(c3+grid%r(1)*c4))))*grid%r(1)**l1
rr2=(1.0_dp+grid%r(2)*(c1+grid%r(2)*(c2+grid%r(2)*(c3+grid%r(2)*c4))))*grid%r(2)**l1
y(1)=rr1/grid%sqr(1)
y(2)=rr2/grid%sqr(2)

return
end subroutine start_scheq



!
! Copyright (C) 2004-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE lschps (mode, z, eps, grid, nin, n, l, e, v, u, nstop)
  !
  ! integrates radial pauli-type scalar-relativistic equation
  ! on a logarithmic grid
  ! modified routine to be used in finding norm-conserving
  ! pseudopotential
  !
  ! on input:
  !   mode = 1 find energy and wavefunction of bound states,
  !            scalar-relativistic (all-electron)
  !   mode = 2 find energy and wavefunction of bound state,
  !            nonrelativistic (pseudopotentials)
  !   mode = 3 fixed-energy calculation, for logarithmic derivatives
  !   mode = 4 find energy which produces a specified logarithmic
  !            derivative (nonrelativistic, pseudopotentials)
  !   mode = 5 is for pseudopotential to produce wavefunction beyond
  !            radius used for pseudopotential construction
  !   z    = atomic number
  !   eps  = convergence factor: eiganvalue is considered converged if
  !          the correction to eigenvalue is smaller in magnitude than
  !          eps times the magnitude of the current guess
  !   grid = structure containing radial grid information
  !   l, n = main and angular quantum numbers
  !   e    = starting estimate of the energy (mode=1,2)
  !          fixed energy at which the wavefctn is calculated (mode=3,4)
  !   v(i) = self-consistent potential
  !   nin  = integration up to r(nin) (mode=3,4,5)
  !
  ! on output:
  !   e    = final energy (mode=1,2)
  !   u(i) = radial wavefunction (defined as the radial part of the wavefct
  !          multiplied by r)
  !   nstop= 0 if regular termination, 1 otherwise
  !   nin  = last grid point for which the wavefct is calculated (mode=1,2)
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : c_au
  USE radial_grids,   ONLY : radial_grid_type
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER, INTENT (in) :: mode, n, l
  real(DP), INTENT(in) :: z, eps
  TYPE (radial_grid_type), INTENT(in) :: grid
  real(DP), INTENT(in) :: v(grid%mesh)
  INTEGER, INTENT(inout) :: nin
  real(DP), INTENT(inout) :: e
  INTEGER, INTENT(out) :: nstop
  real (DP), INTENT(out) :: u(grid%mesh)
  !
  ! local variables
  !
  INTEGER, PARAMETER :: maxter=200
  real(DP), EXTERNAL:: aei, aeo, aii, aio
  ! arrays  used as work space
  real(DP),ALLOCATABLE :: up(:),upp(:),cf(:),dv(:),fr(:),frp(:)
  real(DP):: al, als, cn
  real(DP):: de, emax, emin
  real(DP):: fss, gamma, ro, sc
  real(DP):: sls, sn, uld, uout,  upin, upout
  real(DP):: xkap
  INTEGER:: i, it, mmax, n_it, node, mch, ierr
  !
  !
  nstop=0
  al   = grid%dx
  mmax = grid%mesh

  ALLOCATE(up(mmax), stat=ierr)
  ALLOCATE(upp(mmax), stat=ierr)
  ALLOCATE(cf(mmax), stat=ierr)
  ALLOCATE(dv(mmax), stat=ierr)
  ALLOCATE(fr(mmax), stat=ierr)
  ALLOCATE(frp(mmax), stat=ierr)

  uld=0.0_dp
  !
  !
  IF(mode == 1 .or. mode == 3) THEN
     !     relativistic calculation
     !     fss=(1.0_dp/137.036_dp)**2
     fss=(1.0_dp/c_au)**2
     IF(l == 0) THEN
        gamma=sqrt(1.0_dp-fss*z**2)
     ELSE
        gamma=(l*sqrt(l**2-fss*z**2) + &
             (l+1)*sqrt((l+1)**2-fss*z**2))/(2*l+1)
     ENDIF
  ELSE
     !     non-relativistic calculation
     fss=1.0e-20_dp
     gamma=l+1
  ENDIF
  !
  sls=l*(l+1)
  !
  ! emin, emax = estimated bounds for e
  !
  IF(mode == 1 .or. mode == 2) THEN
     emax=v(mmax)+sls/grid%r(mmax)**2
     emin=0.0_dp
     DO i=1,mmax
        emin=min(emin,v(i)+sls/grid%r(i)**2)
     ENDDO
     IF(e > emax) e=1.25_dp*emax
     IF(e < emin) e=0.75_dp*emin
     IF(e > emax) e=0.5_dp*(emax+emin)
  ELSEIF(mode == 4) THEN
     emax=e + 10.0_dp
     emin=e - 10.0_dp
  ENDIF
  !
  DO i=1,4
     u(i)=0.0_dp
     up(i)=0.0_dp
     upp(i)=0.0_dp
  ENDDO
  als=al**2
  !
  ! calculate dv/dr for darwin correction
  !
  CALL derv (mmax, al, grid%r, v, dv )
  !
  !     starting of loop on energy for bound state
  !
  DO n_it = 1, maxter
     !
     ! coefficient array for u in differential eq.
     DO i=1,mmax
        cf(i)=als*(sls + (v(i)-e)*grid%r(i)**2)
     ENDDO
     !
     ! find classical turning point for matching
     !
     IF(mode == 1 .or. mode == 2) THEN
        DO i=mmax,2,-1
           IF(cf(i-1) <= 0.0_dp .and. cf(i) > 0.0_dp) THEN
              mch=i
              GOTO 40
           ENDIF
        ENDDO
        !PRINT '('' warning: wfc '',2i2,'' no turning point'')', n, l
        e=0.0_dp
        DO i=1,mmax
           u (i)=0.0_dp
        ENDDO
        nstop=1
        GOTO 999
     ELSE
        mch=nin
     ENDIF
40   CONTINUE

     !  relativistic coefficient arrays for u (fr) and up (frp).
     DO i=1,mmax
        fr(i)=als*(grid%r(i)**2)*0.25_dp*(-fss*(v(i)-e)**2 + &
             fss*dv(i)/ (grid%r(i)*(1.0_dp+0.25_dp*fss*(e-v(i)))))
        frp(i)=-al*grid%r(i)*0.25_dp*fss*dv(i)/(1.0_dp+0.25_dp*fss*(e-v(i)))
     ENDDO
     !
     ! start wavefunction with series
     !
     DO i=1,4
        u(i)=grid%r(i)**gamma
        up(i)=al*gamma*grid%r(i)**gamma
        upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
     ENDDO
     !
     ! outward integration using predictor once, corrector
     ! twice
     node=0
     !
     DO i=4,mch-1
        u(i+1)=u(i)+aeo(up,i)
        up(i+1)=up(i)+aeo(upp,i)
        DO it=1,2
           upp(i+1)=(al+frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
           up(i+1)=up(i)+aio(upp,i)
           u(i+1)=u(i)+aio(up,i)
        ENDDO
        IF(u(i+1)*u(i) <= 0.0_dp) node=node+1
     ENDDO
     !
     uout=u(mch)
     upout=up(mch)
     !
     IF(node-n+l+1 == 0 .or. mode == 3 .or. mode == 5) THEN
        !
        IF(mode == 1 .or. mode == 2) THEN
           !
           ! start inward integration at 10*classical turning
           ! point with simple exponential
           nin=mch+2.3_dp/al
           IF(nin+4 > mmax) nin=mmax-4
           xkap=sqrt(sls/grid%r(nin)**2 + 2.0_dp*(v(nin)-e))
           !
           DO i=nin,nin+4
              u(i)=exp(-xkap*(grid%r(i)-grid%r(nin)))
              up(i)=-grid%r(i)*al*xkap*u(i)
              upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
           ENDDO
           !
           ! integrate inward
           !
           DO i=nin,mch+1,-1
              u(i-1)=u(i)+aei(up,i)
              up(i-1)=up(i)+aei(upp,i)
              DO it=1,2
                 upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
                 up(i-1)=up(i)+aii(upp,i)
                 u(i-1)=u(i)+aii(up,i)
              ENDDO
           ENDDO
           !
           ! scale outside wf for continuity
           sc=uout/u(mch)
           !
           DO i=mch,nin
              up(i)=sc*up(i)
              u (i)=sc*u (i)
           ENDDO
           !
           upin=up(mch)
           !
        ELSE
           !
           upin=uld*uout
           !
        ENDIF
        !
        ! perform normalization sum
        !
        ro=grid%r(1)*exp(-0.5_dp*grid%dx)
        sn=ro**(2.0_dp*gamma+1.0_dp)/(2.0_dp*gamma+1.0_dp)
        !
        DO i=1,nin-3
           sn=sn+al*grid%r(i)*u(i)**2
        ENDDO
        !
        sn=sn + al*(23.0_dp*grid%r(nin-2)*u(nin-2)**2 &
             + 28.0_dp*grid%r(nin-1)*u(nin-1)**2 &
             +  9.0_dp*grid%r(nin  )*u(nin  )**2)/24.0_dp
        !
        ! normalize u
        cn=1.0_dp/sqrt(sn)
        uout=cn*uout
        upout=cn*upout
        upin=cn*upin
        !
        DO i=1,nin
           up(i)=cn*up(i)
           u(i)=cn*u(i)
        ENDDO
        DO i=nin+1,mmax
           u(i)=0.0_dp
        ENDDO
        !
        ! exit for fixed-energy calculation
        !
        IF(mode == 3 .or. mode == 5) GOTO 999

        ! perturbation theory for energy shift
        de=uout*(upout-upin)/(al*grid%r(mch))
        !
        ! convergence test and possible exit
        !
        IF ( abs(de) < max(abs(e),0.2_dp)*eps) GOTO 999
        !
        IF(de > 0.0_dp) THEN
           emin=e
        ELSE
           emax=e
        ENDIF
        e=e+de
        IF(e > emax .or. e < emin) e=0.5_dp*(emax+emin)
        !
     ELSEIF(node-n+l+1 < 0) THEN
        ! too few nodes
        emin=e
        e=0.5_dp*(emin+emax)

     ELSE
        ! too many nodes
        emax=e
        e=0.5_dp*(emin+emax)
     ENDIF
  ENDDO

  !PRINT '('' warning: wfc '',2i2,'' not converged'')', n, l
  u=0.0_dp
  nstop=1
  !
  ! deallocate arrays and exit
  !
999 CONTINUE
  DEALLOCATE(frp)
  DEALLOCATE(fr)
  DEALLOCATE(dv)
  DEALLOCATE(cf)
  DEALLOCATE(upp)
  DEALLOCATE(up)
  RETURN

END SUBROUTINE lschps


FUNCTION aei(y,j)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER j
  real(DP):: y(j+3), aei
  !
  aei=-(4.16666666667e-2_dp)*(55.0_dp*y(j)-59.0_dp*y(j+1) &
       +37.0_dp*y(j+2)-9.0_dp*y(j+3))
  RETURN
END FUNCTION aei
!
! adams extrapolation and interpolation formulas for
! outward and inward integration, abramowitz and stegun, p. 896
FUNCTION aeo(y,j)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER:: j
  real(DP):: y(j), aeo
  !
  aeo=(4.16666666667e-2_dp)*(55.0_dp*y(j)-59.0_dp*y(j-1) &
       +37.0_dp*y(j-2)-9.0_dp*y(j-3))
  RETURN
END FUNCTION aeo
!
FUNCTION aii(y,j)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER:: j
  real(DP) :: y(j+2), aii
  !
  aii=-(4.16666666667e-2_dp)*(9.0_dp*y(j-1)+19.0_dp*y(j) &
       -5.0_dp*y(j+1)+y(j+2))
  RETURN
END FUNCTION aii
!
FUNCTION aio(y,j)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: j
  real(DP):: y(j+1), aio
  !
  aio=(4.16666666667e-2_dp)*(9.0_dp*y(j+1)+19.0_dp*y(j) &
       -5.0_dp*y(j-1)+y(j-2))
  RETURN
END FUNCTION aio

SUBROUTINE derV (mmax,al,r,v,dv)
  ! dv = dv/dr
  USE kinds, ONLY : dp
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: mmax
  REAL(dp), INTENT(in) :: al, r(mmax), v(mmax)
  REAL(dp), INTENT(out):: dv(mmax)
  !
  INTEGER :: i
  !
  dv(1)=(-50.0_dp*v(1)+96.0_dp*v(2)-72.0_dp*v(3)+32.0_dp*v(4) &
       -6.0_dp*v(5))/(24.0_dp*al*r(1))
  dv(2)=(-6.0_dp*v(1)-20.0_dp*v(2)+36.0_dp*v(3)-12.0_dp*v(4) &
       +2.0_dp*v(5))/(24.0_dp*al*r(2))
  !
  DO i=3,mmax-2
     dv(i)=(2.0_dp*v(i-2)-16.0_dp*v(i-1)+16.0_dp*v(i+1) &
          -2.0_dp*v(i+2))/(24.0_dp*al*r(i))
  ENDDO
  !
  dv(mmax-1)=( 3.0_dp*v(mmax)+10.0_dp*v(mmax-1)-18.0_dp*v(mmax-2)+ &
       6.0_dp*v(mmax-3)-v(mmax-4))/(12.0_dp*al*r(mmax-1))
  dv(mmax)=( 25.0_dp*v(mmax)-48.0_dp*v(mmax-1)+36.0_dp*v(mmax-2)-&
       16.0_dp*v(mmax-3)+3.0_dp*v(mmax-4))/(12.0_dp*al*r(mmax))
  !
  RETURN
  !
END SUBROUTINE derV


!
! Copyright (C) 2004-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine vxcgc ( ndm, mesh, nspin, r, r2, rho, rhoc, vgc, egc, &
     tau, vtau, iflag)
  !---------------------------------------------------------------
  !
  !
  !     This routine computes the exchange and correlation potential and
  !     energy to be added to the local density, to have the first
  !     gradient correction.
  !     In input the density is rho(r) (multiplied by 4*pi*r2).
  !
  !     The units of the potential are Ry.
  !
  use kinds, only : DP
  use constants, only : fpi, e2
  use funct, only : gcxc, gcx_spin, gcc_spin, dft_is_meta, xc
  implicit none
  integer,  intent(in) :: ndm,mesh,nspin,iflag
  real(DP), intent(in) :: r(mesh), r2(mesh), rho(ndm,2), rhoc(ndm)
  real(DP), intent(out):: vgc(ndm,2), egc(ndm)
  real(DP), intent(in) :: tau(ndm,2)
  real(DP), intent(out):: vtau(mesh)

  integer :: i, is, ierr
  real(DP) :: sx,sc,v1x,v2x,v1c,v2c
  real(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw
  real(DP) :: v3x, v3c, de_cc, dv1_cc,dv2_cc
  real(DP) :: segno, arho
  real(DP) :: rh, zeta, grh2, grho2(2)
  real(DP),parameter :: eps=1.e-12_dp

  real(DP), allocatable :: grho(:,:), h(:,:), dh(:), rhoaux(:,:)
  !
  !      First compute the charge and the charge gradient, assumed  
  !      to have spherical symmetry. The gradient is the derivative of
  !      the charge with respect to the modulus of r. 
  !
  allocate(rhoaux(mesh,2),stat=ierr)
  allocate(grho(mesh,2),stat=ierr)
  allocate(h(mesh,2),stat=ierr)
  allocate(dh(mesh),stat=ierr)

  egc=0.0_dp
  vgc=0.0_dp

  do is=1,nspin
     do i=1, mesh
        rhoaux(i,is)=(rho(i,is)+rhoc(i)/nspin)/fpi/r2(i)
     enddo
     call radial_gradient(rhoaux(1,is),grho(1,is),r,mesh,iflag)
  enddo

  if (nspin.eq.1) then
     !
     IF ( dft_is_meta ()  ) THEN
        !
        !  meta-GGA case
        !
        ! for core correction - not implemented
        de_cc = 0.0_dp
        dv1_cc= 0.0_dp
        dv2_cc= 0.0_dp
        !
	vtau(:) = 0.0_dp
        ! 
       do i=1,mesh
           arho=abs(rhoaux(i,1)) 
           segno=sign(1.0_dp,rhoaux(i,1))
           if (arho.gt.eps.and.abs(grho(i,1)).gt.eps) then
!
! currently there is a single meta-GGA implemented (tpss)
! that calculates all needed terms (LDA, GGA, metaGGA)
!
             call tpsscxc ( arho, grho(i,1)**2, tau(i,1)+tau(i,2), &
                   sx, sc, v1x, v2x, v3x, v1c, v2c, v3c )
              !
              egc(i)=sx+sc+de_cc
              vgc(i,1)= v1x+v1c + dv1_cc
              h(i,1)  =(v2x+v2c)*grho(i,1)*r2(i)
              vtau(i) = v3x+v3c
          else
              vgc(i,1)=0.0_dp
              egc(i)=0.0_dp
              h(i,1)=0.0_dp
              vtau(i)=0.0_dp
           endif
        end do

     ELSE
        !
        !     GGA case
        !
        do i=1,mesh
           arho=abs(rhoaux(i,1)) 
           segno=sign(1.0_dp,rhoaux(i,1))
           if (arho.gt.eps.and.abs(grho(i,1)).gt.eps) then
              call gcxc(arho,grho(i,1)**2,sx,sc,v1x,v2x,v1c,v2c)
              egc(i)=(sx+sc)*segno
              vgc(i,1)= v1x+v1c
              h(i,1)  =(v2x+v2c)*grho(i,1)*r2(i)
              !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
              !                          rho(i,1), grho(i,1)**2,  &
              !                          vgc(i,1),h(i,1)
           else
              vgc(i,1)=0.0_dp
              egc(i)=0.0_dp
              h(i,1)=0.0_dp
           endif
        end do
     END IF
  else
     !
     !   this is the \sigma-GGA case
     !       
     do i=1,mesh
        !
        !  NB: the special or wrong cases where one or two charges 
        !      or gradients are zero or negative must
        !      be detected within the gcxc_spin routine
        !
        !    spin-polarised case
        !
        do is = 1, nspin
           grho2(is)=grho(i,is)**2
        enddo

        call gcx_spin (rhoaux(i, 1), rhoaux(i, 2), grho2(1), grho2(2), &
             sx, v1xup, v1xdw, v2xup, v2xdw)
        rh = rhoaux(i, 1) + rhoaux(i, 2)
        if (rh.gt.eps) then
           zeta = (rhoaux (i, 1) - rhoaux (i, 2) ) / rh
           grh2 = (grho (i, 1) + grho (i, 2) ) **2 
           call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
        else
           sc = 0.0_dp
           v1cup = 0.0_dp
           v1cdw = 0.0_dp
           v2c = 0.0_dp
        endif

        egc(i)=sx+sc
        vgc(i,1)= v1xup+v1cup
        vgc(i,2)= v1xdw+v1cdw
        h(i,1)  =((v2xup+v2c)*grho(i,1)+v2c*grho(i,2))*r2(i)
        h(i,2)  =((v2xdw+v2c)*grho(i,2)+v2c*grho(i,1))*r2(i)
        !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
        !                          rho(i,1)*2.0_dp, grho(i,1)**2*4.0_dp, &
        !                          vgc(i,1),  h(i,2)
     enddo
  endif
  !     
  !     We need the gradient of h to calculate the last part of the exchange
  !     and correlation potential.
  !     
  do is=1,nspin
     call radial_gradient(h(1,is),dh,r,mesh,iflag)
     !
     !     Finally we compute the total exchange and correlation energy and
     !     potential. We put the original values on the charge and multiply
     !     by e^2 = two to have as output Ry units.

     do i=1, mesh
        vgc(i,is)=vgc(i,is)-dh(i)/r2(i)
        vgc(i,is)=e2*vgc(i,is)
        if (is.eq.1) egc(i)=e2*egc(i)
        !            if (is.eq.1.and.i.lt.4) write(6,'(3f20.12)') &
        !                                      vgc(i,1)
     enddo
  enddo
  IF ( dft_is_meta () ) vtau(:) = e2*vtau(:)

  deallocate(dh)
  deallocate(h)
  deallocate(grho)
  deallocate(rhoaux)

  return
end subroutine vxcgc

