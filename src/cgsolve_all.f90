!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define USE_NEW_ROUTINE

#ifdef USE_NEW_ROUTINE
!----------------------------------------------------------------------
SUBROUTINE cgsolve_all(h_psi, cg_psi, e, d0psi, dpsi, h_diag, ndmx, ndim, &
                       ethr, ik, kter, conv_root, anorm, nbnd, npol)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                 ( h - e + Q ) * dpsi = d0psi                      (1)
  !
  !                 where h is a complex hermitean matrix, e is a real sca
  !                 dpsi and d0psi are complex vectors
  !
  !     on input:
  !                 h_psi    EXTERNAL  name of a subroutine:
  !                          h_psi(ndim,psi,psip)
  !                          Calculates  H*psi products.
  !                          Vectors psi and psip should be dimensined
  !                          (ndmx,nvec). nvec=1 is used!
  !
  !                 cg_psi   EXTERNAL  name of a subroutine:
  !                          g_psi(ndmx,ndim,notcnv,psi,e)
  !                          which calculates (h-e)^-1 * psi, with
  !                          some approximation, e.g. (diag(h)-e)
  !
  !                 e        real     unperturbed eigenvalue.
  !
  !                 dpsi     contains an estimate of the solution
  !                          vector.
  !
  !                 d0psi    contains the right hand side vector
  !                          of the system.
  !
  !                 ndmx     integer row dimension of dpsi, ecc.
  !
  !                 ndim     integer actual row dimension of dpsi
  !
  !                 ethr     real     convergence threshold. solution
  !                          improvement is stopped when the error in
  !                          eq (1), defined as l.h.s. - r.h.s., becomes
  !                          less than ethr in norm.
  !
  !     on output:  dpsi     contains the refined estimate of the
  !                          solution vector.
  !
  !                 d0psi    is corrupted on exit
  !
  !   revised 2012-2013 by N. Varini and D. Ceresoli
  !   revised (extensively)       6 Apr 1997 by A. Dal Corso & F. Mauri
  !   revised (to reduce memory) 29 May 2004 by S. de Gironcoli
  !
  USE kinds,          ONLY: dp
  USE mp_pools,       ONLY: intra_pool_comm, me_pool
  USE mp,             ONLY: mp_sum
#ifdef __BANDS
  USE mp_bands,       ONLY: intra_bgrp_comm, inter_bgrp_comm, me_bgrp
  USE gipaw_module,   ONLY: ibnd_start, ibnd_end
#endif
  !-- parameters ------------------------------------------------------
  implicit none
  integer :: ndmx, &      ! input: the maximum dimension of the vectors
             ndim, &      ! input: the actual dimension of the vectors
             kter, &      ! output: counter on iterations
             nbnd, &      ! input: the number of bands
             npol, &      ! input: number of components of the wavefunctions
             ik           ! input: the k point
  real(dp) :: e(nbnd), &  ! input: the actual eigenvalue
              anorm,   &  ! output: the norm of the error in the solution
              h_diag(ndmx*npol,nbnd), & ! input: an estimate of ( H - \epsilon )
              ethr        ! input: the required precision
  complex(dp) :: dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
                 d0psi (ndmx*npol, nbnd)   ! input: the known term
  logical :: conv_root    ! output: if true the root is converged
  external h_psi          ! input: the routine computing h_psi
  external cg_psi         ! input: the routine computing cg_psi

  !-- local variabbles ------------------------------------------------
  integer, parameter :: maxter = 200        ! the maximum number of iterations
  integer :: iter, ibnd                     ! counters on iteration, bands
  integer, allocatable :: conv(:)           ! if 1 the root is converged
  complex(dp), allocatable :: g(:,:), &     ! the gradient of psi
                              t(:,:), &     ! the preconditioned gradient
                              h(:,:), &     ! the delta gradient
                              hold(:,:)     ! the conjugate gradient, work space
  complex(dp) :: dcgamma                    ! the ratio between rho
  complex(dp) :: dclambda                   ! step length
  real(dp), allocatable :: rho(:), rhoold(:), eu(:), a(:), c(:)
  real(dp) :: kter_eff
  complex(dp), external :: zdotc 
#ifndef __BANDS
  integer :: ibnd_start, ibnd_end
  ibnd_start = 1
  ibnd_end = nbnd
#endif

  call start_clock ('cgsolve')

  allocate (g(ndmx*npol,nbnd), t(ndmx*npol,nbnd), h(ndmx*npol,nbnd), hold(ndmx*npol,nbnd))    
  allocate (a(nbnd), c(nbnd))    
  allocate (conv(nbnd))    
  allocate (rho(nbnd),rhoold(nbnd))    
  allocate (eu(nbnd))    

  kter_eff = 0.d0
  conv(1:nbnd) = 0
  g = (0.d0,0.d0)
  t = (0.d0,0.d0)
  h = (0.d0,0.d0)
  hold = (0.d0,0.d0)
  rho = 0.d0

  !--------------------------------------------------------------------
  ! Start CG interations
  !--------------------------------------------------------------------
  do iter = 1, maxter
     ! compute the gradient. can reuse information from previous step
     if (iter == 1) then
#ifdef __BANDS
        call h_psi (ndim, dpsi, g, e, ik, nbnd, ibnd_start, ibnd_end)
#else
        call h_psi (ndim, dpsi, g, e, ik, nbnd)
#endif
        do ibnd = ibnd_start, ibnd_end
           call zaxpy(ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
        enddo

        if (npol == 2) then
           do ibnd = ibnd_start, ibnd_end
              call zaxpy(ndim, (-1.d0,0.d0), d0psi(ndmx+1,ibnd), 1, g(ndmx+1,ibnd), 1)
           enddo
        END IF
     endif

     ! compute preconditioned residual vector and convergence check
     do ibnd = ibnd_start, ibnd_end
        call zcopy (ndmx*npol, g (1, ibnd), 1, h (1, ibnd), 1)
        call cg_psi(ndmx, ndim, 1, h(1,ibnd), h_diag(1,ibnd) )
        rho(ibnd) = zdotc(ndmx*npol, h(1,ibnd), 1, g(1,ibnd), 1)
     enddo
     kter_eff = kter_eff + DBLE (ibnd_end-ibnd_start+1) / DBLE (nbnd)

#ifdef __MPI
#  ifdef __BANDS
     call mp_sum(rho(ibnd_start:ibnd_end), intra_bgrp_comm)
#  else
     call mp_sum(rho(1:nbnd) , intra_pool_comm)
#  endif
#endif

     conv_root = .true.
     anorm = 0.d0
     do ibnd = ibnd_start, ibnd_end
        if (conv(ibnd) == 0) then
           if (sqrt(rho(ibnd)) < ethr) conv(ibnd) = 1
        endif
        conv_root = conv_root .and. (conv(ibnd) == 1)
        anorm = anorm + sqrt(rho(ibnd)) / dble(ibnd_end-ibnd_start+1)
     enddo
     if (conv_root) goto 100

     ! compute the step direction h, conjugate it to previous step
     do ibnd = ibnd_start, ibnd_end
        if (conv(ibnd) == 0) then
           ! change sign to h 
           call dscal(2*ndmx*npol, -1.d0, h(1,ibnd), 1)
           if (iter /= 1) then
              dcgamma = rho(ibnd) / rhoold(ibnd)
              call zaxpy(ndmx*npol, dcgamma, hold(1,ibnd), 1, h(1,ibnd), 1)
           endif
           ! here hold is used as auxiliary vector in order to efficiently compute t = A*h
           ! it is later set to the current (becoming old) value of h 
           call zcopy(ndmx*npol, h(1,ibnd), 1, hold(1,ibnd), 1)
           eu(ibnd) = e(ibnd)
        endif
     enddo

     ! compute t = A*h
#ifdef __BANDS
     call h_psi (ndim, hold, t, eu, ik, nbnd, ibnd_start, ibnd_end)
#else
     call h_psi (ndim, hold, t, eu, ik, nbnd)
#endif

     ! compute the coefficients a and c for the line minimization
     ! compute step length lambda
     do ibnd = ibnd_start, ibnd_end
        if (conv(ibnd) == 0) then
           a(ibnd) = zdotc(ndmx*npol, h(1,ibnd), 1, g(1,ibnd), 1)
           c(ibnd) = zdotc(ndmx*npol, h(1,ibnd), 1, t(1,ibnd), 1)
        endif
     enddo

#ifdef __MPI
#  ifdef __BANDS
     call mp_sum(a(ibnd_start:lbnd), intra_bgrp_comm)
     call mp_sum(c(ibnd_start:lbnd), intra_bgrp_comm)
#  else
     call mp_sum(a(1:nbnd), intra_pool_comm)
     call mp_sum(c(1:nbnd), intra_pool_comm)
#  endif
#endif

     do ibnd = ibnd_start, ibnd_end
        if (conv(ibnd) ==.0) then
           dclambda = cmplx(-a(ibnd)/c(ibnd), 0.d0, kind=dp)
           ! move to new position
           call zaxpy(ndmx*npol, dclambda, h(1,ibnd), 1, dpsi(1,ibnd), 1)
           ! update to get the gradient (g=g+lam)
           call zaxpy(ndmx*npol, dclambda, t(1,ibnd), 1, g(1,ibnd), 1)
           !  save current (now old) h and rho for later use
           call zcopy(ndmx*npol, h(1,ibnd), 1, hold(1,ibnd), 1)
           rhoold(ibnd) = rho(ibnd)
        endif
     enddo
  enddo
  !--------------------------------------------------------------------
  ! End of CG interations
  !--------------------------------------------------------------------
100 continue

#ifdef __BANDS
   do ibnd = ibnd_start, ibnd_end
     if (conv(ibnd) == 0 .and. me_bgrp == 0) &
#else
   do ibnd = 1, nbnd
     if (conv(ibnd) == 0 .and. me_pool == 0) &
#endif
       write(*,'(5x,"ik",i4," ibnd",i4,4x,"cgsolve_all: root not converged")') ik, ibnd
   enddo

#ifdef __BANDS
  call mp_sum(kter_eff, inter_bgrp_comm)
  call mp_sum(anorm, inter_bgrp_comm)
#endif
  kter = kter_eff

  deallocate (eu, rho, rhoold, conv, a, c, g, t, h, hold)

  call stop_clock ('cgsolve')

END SUBROUTINE cgsolve_all


#else
!
! In my opinion, this routine suffers of premature stopping, especially
! for large systems. I don't understand the logic of sorting of the
! converged bands (DC, Oct 2013)
!
!----------------------------------------------------------------------
subroutine cgsolve_all (h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                 ( h - e + Q ) * dpsi = d0psi                      (1)
  !
  !                 where h is a complex hermitean matrix, e is a real sca
  !                 dpsi and d0psi are complex vectors
  !
  !     on input:
  !                 h_psi    EXTERNAL  name of a subroutine:
  !                          h_psi(ndim,psi,psip)
  !                          Calculates  H*psi products.
  !                          Vectors psi and psip should be dimensined
  !                          (ndmx,nvec). nvec=1 is used!
  !
  !                 cg_psi   EXTERNAL  name of a subroutine:
  !                          g_psi(ndmx,ndim,notcnv,psi,e)
  !                          which calculates (h-e)^-1 * psi, with
  !                          some approximation, e.g. (diag(h)-e)
  !
  !                 e        real     unperturbed eigenvalue.
  !
  !                 dpsi     contains an estimate of the solution
  !                          vector.
  !
  !                 d0psi    contains the right hand side vector
  !                          of the system.
  !
  !                 ndmx     integer row dimension of dpsi, ecc.
  !
  !                 ndim     integer actual row dimension of dpsi
  !
  !                 ethr     real     convergence threshold. solution
  !                          improvement is stopped when the error in
  !                          eq (1), defined as l.h.s. - r.h.s., becomes
  !                          less than ethr in norm.
  !
  !     on output:  dpsi     contains the refined estimate of the
  !                          solution vector.
  !
  !                 d0psi    is corrupted on exit
  !
  !   revised (extensively)       6 Apr 1997 by A. Dal Corso & F. Mauri
  !   revised (to reduce memory) 29 May 2004 by S. de Gironcoli
  !
  USE kinds, only : DP
  USE mp_pools,  ONLY: intra_pool_comm, me_pool
  USE mp,        ONLY: mp_sum
#ifdef __BANDS
  USE mp_bands,  ONLY: intra_bgrp_comm, inter_bgrp_comm, me_bgrp
  USE gipaw_module, ONLY: ibnd_start, ibnd_end
#endif

  implicit none
  !
  !   first the I/O variables
  !
  integer :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             npol, & ! input: number of components of the wavefunctions
             ik      ! input: the k point

  real(DP) :: &
             e(nbnd), & ! input: the actual eigenvalue
             anorm,   & ! output: the norm of the error in the solution
             h_diag(ndmx*npol,nbnd), & ! input: an estimate of ( H - \epsilon )
             ethr       ! input: the required precision

  complex(DP) :: &
             dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx*npol, nbnd)   ! input: the known term

  logical :: conv_root ! output: if true the root is converged
  external h_psi       ! input: the routine computing h_psi
  external cg_psi      ! input: the routine computing cg_psi
  !
  !  here the local variables
  !
  integer, parameter :: maxter = 200
  ! the maximum number of iterations
  integer :: iter, ibnd, lbnd
  ! counters on iteration, bands
  integer , allocatable :: conv (:)
  ! if 1 the root is converged

  complex(DP), allocatable :: g (:,:), t (:,:), h (:,:), hold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  complex(DP) ::  dcgamma, dclambda
  !  the ratio between rho
  !  step length
  complex(DP), external :: zdotc 
  !  the scalar product
  real(DP), allocatable :: rho (:), rhoold (:), eu (:), a(:), c(:)
  ! the residue
  ! auxiliary for h_diag
  real(DP) :: kter_eff
  ! account the number of iterations with b
  ! coefficient of quadratic form
  !
  call start_clock ('cgsolve')
  allocate ( g(ndmx*npol,nbnd), t(ndmx*npol,nbnd), h(ndmx*npol,nbnd), &
             hold(ndmx*npol ,nbnd) )    
  allocate (a(nbnd), c(nbnd))    
  allocate (conv ( nbnd))    
  allocate (rho(nbnd),rhoold(nbnd))    
  allocate (eu (  nbnd))    
  !      WRITE( stdout,*) g,t,h,hold

  kter_eff = 0.d0
  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo
  g=(0.d0,0.d0)
  t=(0.d0,0.d0)
  h=(0.d0,0.d0)
  hold=(0.d0,0.d0)
  do iter = 1, maxter
     !
     !    compute the gradient. can reuse information from previous step
     !
     if (iter == 1) then
#ifdef __BANDS
        call h_psi (ndim, dpsi, g, e, ik, nbnd, ibnd_start, ibnd_end)
        do ibnd = ibnd_start, ibnd_end
#else
        call h_psi (ndim, dpsi, g, e, ik, nbnd)
        do ibnd = 1, nbnd
#endif
           call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
        enddo
        IF (npol==2) THEN
#ifdef __BANDS
           do ibnd = ibnd_start, ibnd_end
#else
           do ibnd = 1, nbnd
#endif
              call zaxpy (ndim, (-1.d0,0.d0), d0psi(ndmx+1,ibnd), 1, &
                                              g(ndmx+1,ibnd), 1)
           enddo
        END IF
     endif
     !
     !    compute preconditioned residual vector and convergence check
     !
#ifdef __BANDS
     lbnd = ibnd_start - 1
     do ibnd = ibnd_start, ibnd_end
#else
     lbnd = 0
     do ibnd = 1, nbnd
#endif
        if (conv (ibnd) .eq.0) then
           lbnd = lbnd+1
           call zcopy (ndmx*npol, g (1, ibnd), 1, h (1, ibnd), 1)
           call cg_psi(ndmx, ndim, 1, h(1,ibnd), h_diag(1,ibnd) )
           rho(lbnd) = zdotc (ndmx*npol, h(1,ibnd), 1, g(1,ibnd), 1)
        endif
     enddo

#ifdef __BANDS
     kter_eff = kter_eff + DBLE (lbnd-ibnd_start+1) / DBLE (nbnd)
#else
     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
#endif

#ifdef __MPI
#  ifdef __BANDS
     call mp_sum(  rho(ibnd_start:lbnd), intra_bgrp_comm )
#  else
     call mp_sum(  rho(1:lbnd) , intra_pool_comm )
#  endif
#endif

#ifdef __BANDS
     do ibnd = ibnd_end, ibnd_start, -1
#else
     do ibnd = nbnd, 1, -1
#endif
        if (conv(ibnd).eq.0) then
           rho(ibnd)=rho(lbnd)
           lbnd = lbnd - 1
           anorm = sqrt(rho(ibnd))
           if (anorm.lt.ethr) conv (ibnd) = 1
        endif
     enddo

     conv_root = .true.
     anorm = 0.d0
#ifdef __BANDS
     do ibnd = ibnd_start, ibnd_end
#else
     do ibnd = 1, nbnd
#endif
        conv_root = conv_root .and. (conv(ibnd) .eq.1)
        anorm = anorm + sqrt(rho(ibnd)) / dble(nbnd)
     enddo

     if (conv_root) goto 100
     !
     !        compute the step direction h. Conjugate it to previous step
     !
#ifdef __BANDS
     lbnd = ibnd_start - 1 
     do ibnd = ibnd_start, ibnd_end
#else
     lbnd = 0
     do ibnd = 1, nbnd
#endif
        if (conv (ibnd) .eq.0) then
!
!          change sign to h 
!
           call dscal (2 * ndmx * npol, - 1.d0, h (1, ibnd), 1)
           if (iter.ne.1) then
              dcgamma = rho (ibnd) / rhoold (ibnd)
              call zaxpy (ndmx*npol, dcgamma, hold (1, ibnd), 1, h (1, ibnd), 1)
           endif

!
! here hold is used as auxiliary vector in order to efficiently compute t = A*h
! it is later set to the current (becoming old) value of h 
!
           lbnd = lbnd+1
           call zcopy (ndmx*npol, h (1, ibnd), 1, hold (1, lbnd), 1)
           eu (lbnd) = e (ibnd)
        endif
     enddo
     !
     !        compute t = A*h
     !
#ifdef __BANDS
     call h_psi (ndim, hold, t, eu, ik, nbnd, ibnd_start, lbnd)
#else
     call h_psi (ndim, hold, t, eu, ik, lbnd)
#endif
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
#ifdef __BANDS
     lbnd = ibnd_start - 1
     do ibnd = ibnd_start, ibnd_end
#else
     lbnd=0
     do ibnd = 1, nbnd
#endif
        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1
           a(lbnd) = zdotc (ndmx*npol, h(1,ibnd), 1, g(1,ibnd), 1)
           c(lbnd) = zdotc (ndmx*npol, h(1,ibnd), 1, t(1,lbnd), 1)
        end if
     end do
#ifdef __MPI
#  ifdef __BANDS
     call mp_sum(  a(ibnd_start:lbnd), intra_bgrp_comm )
     call mp_sum(  c(ibnd_start:lbnd), intra_bgrp_comm )
#  else
     call mp_sum(  a(1:lbnd), intra_pool_comm )
     call mp_sum(  c(1:lbnd), intra_pool_comm )
#  endif
#endif

#ifdef __BANDS
     lbnd = ibnd_start - 1 
     do ibnd = ibnd_start, ibnd_end
#else
     lbnd=0
     do ibnd = 1, nbnd
#endif
        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1
           dclambda = CMPLX( - a(lbnd) / c(lbnd), 0.d0,kind=DP)
           !
           !    move to new position
           !
           call zaxpy (ndmx*npol, dclambda, h(1,ibnd), 1, dpsi(1,ibnd), 1)
           !
           !    update to get the gradient
           !
           !g=g+lam
           call zaxpy (ndmx*npol, dclambda, t(1,lbnd), 1, g(1,ibnd), 1)
           !
           !    save current (now old) h and rho for later use
           ! 
           call zcopy (ndmx*npol, h(1,ibnd), 1, hold(1,ibnd), 1)
           rhoold (ibnd) = rho (ibnd)
        endif
     enddo
  enddo
100 continue

#ifdef __BANDS
   do ibnd = ibnd_start, ibnd_end
     if (conv(ibnd) == 0 .and. me_bgrp == 0) &
#else
   do ibnd = 1, nbnd
     if (conv(ibnd) == 0 .and. me_pool == 0) &
#endif
       write(*,'(5x,"ik",i4," ibnd",i4,4x,"cgsolve_all: root not converged")') ik, ibnd
   enddo

#ifdef __BANDS
  call mp_sum(kter_eff, inter_bgrp_comm)
  call mp_sum(anorm, inter_bgrp_comm)
#endif
  kter = kter_eff
  deallocate (eu)
  deallocate (rho, rhoold)
  deallocate (conv)
  deallocate (a,c)
  deallocate (g, t, h, hold)

  call stop_clock ('cgsolve')
  return
end subroutine cgsolve_all
#endif
