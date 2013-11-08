!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#ifdef __BANDS 
!-----------------------------------------------------------------------
SUBROUTINE init_parallel_over_band(comm, nbnd)
  !-----------------------------------------------------------------------
  !
  ! ... Setup band indexes for band parallelization
  !
  USE gipaw_module, ONLY : ibnd_start, ibnd_end
  USE io_global,    ONLY : stdout
  USE mp_bands,     ONLY : nbgrp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: comm, nbnd 

  INTEGER :: mp_size, mp_rank, ierror, rest, k

  call mpi_comm_rank(comm, mp_rank, ierror)
  call mpi_comm_size(comm, mp_size, ierror)
  rest = mod(nbnd, mp_size)
  k = int(nbnd/mp_size)
    
  if (k.ge.1) then
     if (rest > mp_rank) then
        ibnd_start = (mp_rank)*k + (mp_rank+1)
        ibnd_end  =  (mp_rank+1)*k + (mp_rank+1)
     else
        ibnd_start = (mp_rank)*k + rest + 1
        ibnd_end  =  (mp_rank+1)*k + rest
     endif
  else
     ibnd_start = 1
     ibnd_end = nbnd
  endif   

  write(stdout,'(5X,''============================================================'')')
  write(stdout,'(5X,''  EXPERIMENTAL BAND PARALLELIZATION OVER '',I2,'' GROUPS'')') nbgrp
  write(stdout,'(5X,''============================================================'')')
  write(stdout,*)
END SUBROUTINE init_parallel_over_band
#endif


#ifdef __BANDS
!-----------------------------------------------------------------------
SUBROUTINE calbec_bands (npwx, npw, nkb, beta, psi, betapsi, nbnd, ibnd_start, ibnd_end)
  !-----------------------------------------------------------------------
  !
  ! ... matrix times matrix with summation index (k=1,npw) running on
  ! ... G-vectors or PWs : betapsi(i,j) = \sum_k beta^*(i,k) psi(k,j)
  !
  USE kinds,                 only : dp
  USE mp_bands,              only : intra_bgrp_comm
  USE mp,                    only : mp_sum
  IMPLICIT NONE
  integer, intent(in) :: npwx, npw, nkb
  complex(dp), intent(in) :: beta(npwx,nkb), psi(npwx,nbnd)
  complex(dp), intent(out) :: betapsi(nkb,nbnd)
  integer, intent(in) :: nbnd, ibnd_start, ibnd_end

  if (nkb == 0) return

  call start_clock( 'calbec' )

  call ZGEMM('C', 'N', nkb, ibnd_end-ibnd_start+1, npw, (1.d0,0.d0), &
             beta, npwx, psi(1,ibnd_start), npwx, (0.d0,0.d0), betapsi(1,ibnd_start), nkb)

  call mp_sum(betapsi(:,1:nbnd), intra_bgrp_comm)

  call stop_clock( 'calbec' )

  return
END SUBROUTINE calbec_bands
#endif


#ifdef __BANDS
!-----------------------------------------------------------------------
SUBROUTINE s_psi_bands (lda, n, m, psi, spsi, ibnd_start, ibnd_end)
  !-----------------------------------------------------------------------
  USE kinds,                 only : dp
  USE becmod,                only : becp
  USE uspp,                  only : qq, vkb, nkb, okvan
  USE uspp_param,            only : upf, nhm, nh
  USE ions_base,             only : nat, ityp, ntyp => nsp
  IMPLICIT NONE
  ! -- parameters --------------------------------------------------------
  integer, intent(in) :: lda, n, m
  integer, intent(in) :: ibnd_start, ibnd_end
  complex(dp), intent(in) :: psi(lda,m)
  complex(dp), intent(out) :: spsi(lda,m)
  ! -- local variables ---------------------------------------------------
  integer :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
  complex(dp), allocatable :: ps(:,:)

  ! apply identity
  spsi = (0.d0,0.d0) 
  spsi(1:lda, ibnd_start:ibnd_end) = psi(1:lda,ibnd_start:ibnd_end)

  if (nkb == 0 .or. .not. okvan) return

  allocate( ps(nkb,m) )
  ps(:,:) = (0.d0,0.d0) 

  ijkb0 = 0
  do nt = 1, ntyp
    if ( upf(nt)%tvanp ) then
      do na = 1, nat
        if ( ityp(na) == nt ) then

          do ibnd = ibnd_start, ibnd_end
            do jh = 1, nh(nt)
              jkb = ijkb0 + jh
              do ih = 1, nh(nt)
                ikb = ijkb0 + ih
                ps(ikb,ibnd) = ps(ikb,ibnd) + qq(ih,jh,nt) * becp%k(jkb,ibnd)
              end do  ! ikb
            end do  ! jkb
          end do  ! ibnd

          ijkb0 = ijkb0 + nh(nt)
        end if
      end do ! na
 
     else
       do na = 1, nat
         if ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
       end do
    end if
  end do ! nt

  call ZGEMM('N', 'N', n, ibnd_end-ibnd_start+1, nkb, (1.d0,0.d0), vkb, &
             lda, ps(1,ibnd_start), nkb, (1.d0,0.d0), spsi(1,ibnd_start), lda)

  deallocate( ps )
  return
END SUBROUTINE s_psi_bands
#endif


#ifdef __BANDS
!-----------------------------------------------------------------------
subroutine add_vuspsi_bands (lda, n, m, psi, ibnd_start, ibnd_end)
  !-----------------------------------------------------------------------
  USE kinds,                 only : dp
  USE becmod,                only : becp
  USE uspp,                  only : qq, vkb, nkb, deeq
  USE uspp_param,            only : upf, nhm, nh
  USE ions_base,             only : nat, ityp, ntyp => nsp
  USE lsda_mod,              only : current_spin
  IMPLICIT NONE
  ! -- parameters --------------------------------------------------------
  integer, intent(in) :: lda, n, m
  integer, intent(in) :: ibnd_start, ibnd_end
  complex(dp), intent(inout) :: psi(lda,m)
  ! -- local variables ---------------------------------------------------
  integer :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
  complex(dp), allocatable :: ps(:,:)
  
  if ( nkb == 0 ) return

  allocate( ps(nkb,m) )
  ps(:,:) = (0.d0,0.d0) 

  ijkb0 = 0
  do nt = 1, ntyp
    do na = 1, nat
      if ( ityp(na) == nt ) then

        do ibnd = ibnd_start, ibnd_end
          do jh = 1, nh(nt)
            jkb = ijkb0 + jh
            do ih = 1, nh(nt)
              ikb = ijkb0 + ih
                ps(ikb,ibnd) = ps(ikb,ibnd) + deeq(ih,jh,na,current_spin) * becp%k(jkb,ibnd)
            end do  ! ikb
          end do  ! jkb
        end do  ! ibnd

        ijkb0 = ijkb0 + nh(nt)
      end if
    end do  ! na
  end do  ! nt 

  call ZGEMM('N', 'N', n, ibnd_end-ibnd_start+1, nkb, (1.d0,0.d0) , vkb, &
             lda, ps(1,ibnd_start), nkb, (1.d0,0.d0), psi(1,ibnd_start), lda)

  deallocate (ps)
  return
END SUBROUTINE add_vuspsi_bands
#endif

