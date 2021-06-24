! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine gets estimates of PBL top, using 1 of 4 methods:     !
!               method = 0: min of wt flux, buoyancy cases                     !
!               method = 1: uw flux less than critical value, ekman cases      !
!               method = 2: running t average exceeds criterion, neutral cases !
!                           with capping inversions                            !
!               method = 3: maximum gradient in temperature field, finding     !
!                           local zi (Jas et al.)                              !
!                                                                              !
!         Note that if the chosen method uses average statistics, only root    !
!         process need find zi.                                                !
! ============================================================================ !
!
SUBROUTINE pbltop(itop)
!
    USE pars
    USE fields
    USE con_data
    USE con_stats
!
    INCLUDE 'mpif.h'
!
    REAL :: trun(maxnz), gradloc(2,nnx,nny), gradMAX(2,nnx,nny)
!
    EXTERNAL get_zi
!
! --------------------------------------------------------------------------- !
!
  IF (method .LE. 2 .AND. l_root) THEN
    sgn = 1.0
!
    IF (method .LE. 0 .OR. method .GT. 2) THEN
      itop  = 1
      wttot = wtle(1,1)+wtsb(1,1)
      wtmin = wttot*sgn
!
      DO iz=2,nnz
        wttot = (wtle(iz,1)+wtsb(iz,1))*sgn
        IF (wttot .LE. wtmin) THEN
          itop  = iz
          wtmin = wttot
        END IF
      END DO
!
      zi = z(itop)
!
    ELSE IF (method .EQ. 1) THEN
!
      itop = 1
      crit = 0.05
      uwsf = utau*utau
!
      DO iz=1,nnzm1
        uwtot = (uwle(iz) + uwsb(iz))**2 + (vwle(iz) + vwsb(iz))**2
        uwtot = SQRT(uwtot)
        IF (uwtot/uwsf .GT. crit) THEN
          itop=iz
        END IF
      END DO
!
      zi = z(itop)
!
    ELSE IF (method .EQ. 2) THEN
!
      trun(1) = txym(1,1)
!
      DO iz=2,nnz
        weight   = z(iz-1)/z(iz)
        trun(iz) = trun(iz-1)*weight + (1.0-weight)*txym(iz,1)
      END DO
!
      itop  = 1
      tcrit = 0.1
!
      DO iz=2,nnz
        IF (txym(iz,1) .GT. (trun(iz) + tcrit)) THEN
          itop = iz
          GO TO 320
        END IF
      END DO
!
!
 320  CONTINUE
    zi = z(itop)
    END IF
!
    DO iy=1,nny
      DO ix=1,nnx
        gradMAX(2,ix,iy) = zi
      END DO
    END DO
!
! USE GRADIENT METHOD
  ELSE IF (method .EQ. 3) THEN
!
!   SIMILAR TO ZEROING THE STAT ARRAY
    DO iy=1,nny
      DO ix=1,nnx
        gradloc(1,ix,iy) = 0.0
        gradloc(2,ix,iy) = z(iz_min)
      END DO
    END DO
!
!   ACCOUNT FOR ALL Z
    IF (iz_min .LE. ize) THEN
      DO iz = MAX(izs,iz_min),ize
        izp1 = iz + 1
        DO iy=iys,iye
          DO ix=1,nnx
            grad = (t(ix,iy,1,izp1) - t(ix,iy,1,iz))*dzu_i(izp1)
            IF (grad .GT. gradloc(1,ix,iy)) THEN
              gradloc(1,ix,iy) = grad
              gradloc(2,ix,iy) = z(iz)
            END IF
          END DO
        END DO
      END DO
    END IF
!
!   ALTERNATE VERSION DEFINED IN MPI - PASSES 2 REAL8 VARIABLES
    CALL mpi_reduce(gradloc,gradmax,nnx*nny,mpi_2double_precision,  &
                        mpi_maxloc,i_root,mpi_comm_world,ierror)
!
!   FIND AVERAGE ON ROOT PROCESS
    IF (l_root) THEN
      zi_avg = 0.0
      DO iy=1,nny
        DO ix=1,nnx
          zi_avg = zi_avg + gradMAX(2,ix,iy)
        END DO
      END DO
      zi = zi_avg*fnxy
    END IF
  END IF
!
! SEND ZI AVERAGE
  CALL mpi_bcast(zi,1,mpi_real8,i_root,mpi_comm_world,ierr)
  DO iz=1,nnz
    IF (zi .LE. z(iz) .AND. zi .GT. z(iz+1)) itop = iz
  END DO
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
7001 FORMAT (' 7001 in pbltop myid = ',i4,' zi = ',e15.6,' itop = ',i3)
!------------------------------------------------------------------------------!
!
END SUBROUTINE pbltop
