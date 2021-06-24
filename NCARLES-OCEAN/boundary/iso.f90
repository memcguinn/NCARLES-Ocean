! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine solves for the isotropy factor and scales it to      !
!         match for a corresponding height.                                    !
!                                                                              !
!         Lower and upper boundary conditions are used.                        !
! ============================================================================ !
!
SUBROUTINE iso(it)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    REAL :: sfk(1:nnz)
!
! --------------------------------------------------------------------------- !
!
  DO iz=1,nnz
    dfac(iz) = 0.0
    sfk(iz)  = 0.0
  END DO
  DO iz=izs,ize
    dfac(iz) = 1.0
  END DO
!
! SET MATCH EQUAL TO FRACTION OF INITIAL ZI
  nmatch = nnz
  DO i=0,numprocs-1,ncpu_s
    IF (nmatch .GE. iz_s(i) .AND. nmatch .LE. iz_e(i)) myid_newvis = i
  END DO
  DO iz=izs,MIN(ize,nmatch)
    izp1 = iz + 1
    izm1 = iz - 1
    weit = dzw(iz)/(dzw(iz) + dzw(izp1))
    weit1 = 1.0 - weit
!
!   SOLVE FOR FLUCTUATING STRAIN
    DO j=iys,iye
      DO i=1,nnx
        s11 = weit1*ux(i,j,iz)**2 + weit*ux(i,j,izp1)**2
        s22 = weit1*vy(i,j,iz)**2 + weit*vy(i,j,izp1)**2
        wz  = (w(i,j,iz)-w(i,j,izm1))*dzw_i(iz)
        wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
        s33 = weit*wzp**2 + weit1*wz**2
        s12 = weit1*(uy(i,j,iz) + vx(i,j,iz))**2 +         &
              weit*(uy(i,j,izp1) + vx(i,j,izp1))**2
        s13 = (((u(i,j,izp1) - u(i,j,iz) +                 &
              u_mn(iz) - u_mn(izp1))*dzu_i(izp1) +         &
              wx(i,j,iz)))**2
        s23 = (((v(i,j,izp1) - v(i,j,iz) +                 &
              v_mn(iz) - v_mn(izp1))*dzu_i(izp1) +         &
              wy(i,j,iz)))**2
        sfk(iz) = sfk(iz) + 2.0*(s11 + s22 + s33) + s12 + s13 + s23
      END DO
    END DO
    sfk(iz) = sfk(iz)*fnxy
  END DO
!
  CALL mpi_sum_z(sfk,i_root,myid,nnz,1)
  DO iz=izs,MIN(ize,nmatch)
    izp1 = iz + 1
    izm1 = iz - 1
    sfk(iz) = SQRT(sfk(iz))
    smk = SQRT((u_mn(izp1)-u_mn(iz))**2 +  &
          (v_mn(izp1)-v_mn(iz))**2)*ABS(dzu_i(izp1))
    IF (sfk(iz) .LE. 0. .AND. smk .LE. 0.) THEN
      dfac(iz) = 1.0
    ELSE
      dfac(iz) = sfk(iz)/(sfk(iz) + smk)
    END IF
  END DO
!
! RESCALE RATIO TO MATCH HEIGHT AND NESTED GRID, MATCH VALUE AT UPPER BOUNDARY
! OF COURSER GRID
  IF (myid .EQ. myid_newvis) THEN
      dfacm = dfac(nmatch)
  END IF
  CALL mpi_bcast(dfacm,1,mpi_real8,myid_newvis,mpi_comm_world,ierr)
  DO iz=izs,MIN(ize,nmatch)
    dfac(iz) = dfac(iz)/dfacm
    dfac(iz) = AMAX1(dfac(iz), 0.1)
    dfac(iz) = AMIN1(dfac(iz), 1.0)
  END DO
!
! GATHER DFAC ON ALL PROCESSES FOR PRINTING AND USE IN TKE_VIS
! REDUCE AND DIVIDE BY NUMBER OF CPU SLABS
  CALL mpi_sum_z(dfac,i_root,myid,nnz,1)
  fncpu_s = 1.0/FLOAT(ncpu_s)
  DO iz=1,nnz
    dfac(iz) = dfac(iz)*fncpu_s
  END DO
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
6001 FORMAT(' iz = ',i3,' sfk = ',e15.6,' smk = ',e15.6,' dfac = ',e15.6)
6000 FORMAT(' in sr. iso, nmatch = ',i3,/,'iz',5x,'dfac',/,(i3,1x,e15.6))
3001 FORMAT(' iz ',5x,' dfac ',/,(i5,e15.6))
!------------------------------------------------------------------------------!
!
END SUBROUTINE iso
