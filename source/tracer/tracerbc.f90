! ============================================================================ !
! ABOUT:                                                                       !
!         This module uses Zeebe et al. (2001) carbonate chemistry
!
!         Need equispaced z grid cells for *toi* functions to work.
! ============================================================================ !
!
MODULE tracerbc
!
  USE con_data, ONLY: dz,dx,dy
  USE con_stats, ONLY: z, zz
  USE pars, ONLY: iys,iye,izs,ize,izi,nnx,nnz,nny,nscl
  USE fields, ONLY: t
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: flg_debug = 0         ! Write debug file (0/1)
!
  REAL, DIMENSION(nscl) :: tau,             & ! Reaction timescale
                           airval             ! Value of air tracer (only need with asflux on)
!
  INTEGER, DIMENSION(nscl) :: ictype,       & ! Initial condition
                              rmodel,       & ! Reaction model type (0 - no reaction, 3 - carbonate chemistry)
                              rdorg,        & ! Decaying (0), growing (1) reaction
                              rpartner,     & ! Reaction partner (array of iscl, 0 for no coupled tracer)
                              asflux          ! Air-sea flux boundary condition in place (0/1)
!
  INTEGER, DIMENSION(2,nscl) :: bnd
  INTEGER, DIMENSION(2) :: bnds
  INTEGER, DIMENSION(3,nscl) :: point
  INTEGER, DIMENSION(3) :: points
  REAL, DIMENSION(nscl) :: val                 ! Value of initial finite or source point
!
CONTAINS
!
! --------------------------------------------------------------------------- !
!
    SUBROUTINE applytracerbc(it)
!
      INTEGER, INTENT(in) :: it
      INTEGER :: iscl,                        & ! Scalar number (theta => iscl=1)
                 np,                          & ! Width of initial or finite source
                 zt                             ! Upper/left most level
      REAL :: ta,                             &
              vals
!
! ACTIVE TRACERS, TEMPERATURE
      iscl = 1;                                 !
      ictype(iscl) = 0;                         ! Iniital condition
      val(iscl) = 298.15;
      tau(iscl)    = 0;
      asflux(iscl) = 0;
      airval(iscl) = 0;
      np = 0;
      zt = 0;
      rmodel(iscl) = 0;
      bnd(:,iscl) = znptobnd(zt,np);
!
! PASSIVE TRACERS
!   CARBON DIOXIDE
      iscl = 2;
      ictype(iscl) = 1;
      val(iscl) = 7.56903;
      tau(iscl) = 1;
      asflux(iscl) = 0;
      airval(iscl) = 7.56903 * 1.1;             ! 10% Flux in from atmosphere
      np = nnz+2;
      zt = 0;
      rmodel(iscl) = 3;
      bnd(:,iscl) = znptobnd(zt,np);
!   BICARBONATE
      iscl = 3;
      ictype(iscl) = 1;
      val(iscl) = 1.67006e03;
      tau(iscl) = 1;
      asflux(iscl) = 0;
      airval(iscl) = 0;
      np = nnz+2;
      zt = 0;
      rmodel(iscl) = 3;
      bnd(:,iscl) = znptobnd(zt,np);
!   CARBONATE
      iscl = 4;
      ictype(iscl) = 1;
      val(iscl) = 3.14655e02;
      tau(iscl) = 1;
      asflux(iscl) = 0;
      airval(iscl) = 0;
      np = nnz+2;
      zt = 0;
      rmodel(iscl) = 3;
      bnd(:,iscl) = znptobnd(zt,np);
!   BORIC ACID
      iscl = 5;
      ictype(iscl) = 1;
      val(iscl) = 2.96936e02;
      tau(iscl) = 1;
      asflux(iscl) = 0;
      airval(iscl) = 0;
      np = nnz+2;
      zt = 0;
      rmodel(iscl) = 3;
      bnd(:,iscl) = znptobnd(zt,np);
!   TETRAHYDROXYBORATE
      iscl = 6;
      ictype(iscl) = 1;
      val(iscl) = 1.18909e02;
      tau(iscl) = 1;
      asflux(iscl) = 0;
      airval(iscl) = 0;
      np = nnz+2;
      zt = 0;
      rmodel(iscl) = 3;
      bnd(:,iscl) = znptobnd(zt,np);
!   HYDROGEN ION
      iscl = 7;
      ictype(iscl) = 1;
      val(iscl) = 6.30928e-03;
      tau(iscl) = 1;
      asflux(iscl) = 0;
      airval(iscl) = 0;
      np = nnz+2;
      zt = 0;
      rmodel(iscl) = 3;
      bnd(:,iscl) = znptobnd(zt,np);
!   HYDROXIDE
      iscl = 8;
      ictype(iscl) = 1;
      val(iscl) = 9.60492;
      tau(iscl) = 1;
      asflux(iscl) = 0;
      airval(iscl) = 0;
      np = nnz+2;
      zt = 0;
      rmodel(iscl) = 3;
      bnd(:,iscl) = znptobnd(zt,np);
!
!
      DO iscl = 2,nscl
         bnds=bnd(:,iscl); vals=val(iscl); points=point(:,iscl);
         IF (it.EQ.1) THEN
            IF (ictype(iscl).EQ.1) CALL hbndsource(iscl,bnds,vals);
            IF (ictype(iscl).EQ.5) CALL vgradsource(iscl,bnds,vals);
         ENDIF
         bnds = 0;
         vals = 0;
         points = 0;
      ENDDO
!
      IF(flg_debug == 1) THEN
          OPEN(13, file='tracerbc.txt',access='append')
          WRITE(13,'(A)') '------------------------'
          WRITE(13,'(A,i3)') 'RUNNING FOR IT= ',it
          WRITE(13,'(A,f9.6)') 'Z for 5m above H', z(izi)+5.0
          WRITE(13,'(A,f9.6)') 'Z for 5m below H', z(izi)-5.0
          CLOSE(13)
      END IF
    END SUBROUTINE
!
!
! --------------------------------------------------------------------------- !
!
    SUBROUTINE hbndsource(iscl, bnds, vals)
!
      INTEGER, INTENT(in) :: iscl
      INTEGER, INTENT(in), DIMENSION(2) :: bnds
      REAL, INTENT(in) :: vals
      INTEGER :: ix,iy,iz
!
      DO iz=bnds(1),bnds(2)
         DO iy=iys,iye
            DO ix=1,nnx
               IF ((iz >= izs) .AND. (iz <= ize)) THEN
                     t(ix,iy,iscl,iz) = vals
               ENDIF
            END DO
         END DO
      END DO
!
    END SUBROUTINE
!
!
! --------------------------------------------------------------------------- !
!
    SUBROUTINE vgradsource(iscl, bnds, vals)
!
      INTEGER, INTENT(in) :: iscl
      INTEGER, INTENT(in), DIMENSION(2) :: bnds
      REAL, INTENT(in) :: vals
      INTEGER :: ix,iy,iz,zi
!
      zi  = z(bnds(2))
      DO iy=iys,iye
         DO iz=bnds(1),bnds(2)
            DO ix=1,nnx
               IF ((iz >= izs) .AND. (iz <= ize)) THEN
                  t(ix,iy,iscl,iz) = (vals/zi)*(zi-zz(iz))
               ENDIF
            END DO
         END DO
      END DO
!
    END SUBROUTINE
!
!
! --------------------------------------------------------------------------- !
!
    FUNCTION znptobnd(zt,np)
!
      INTEGER, INTENT(in) :: zt
      INTEGER, INTENT(in) :: np
      INTEGER, DIMENSION(2) :: znptobnd
      INTEGER :: iz
!
!     CHECK: DIMENSIONS WHEN SETTING FIRST BOUND
      iz = ztoiz(zt)
      znptobnd(1) = iz - INT((np-1)/2)
!
      IF (znptobnd(1) < 0) THEN
        znptobnd(1) = 0
      END IF
!
!     SET SECOND BOUND BY FIRST
      znptobnd(2) = znptobnd(1) + np -1
!
    END FUNCTION
!
!
! --------------------------------------------------------------------------- !
!
    FUNCTION xnptobnd(xt,np)
!
      INTEGER, INTENT(in) :: xt
      INTEGER, INTENT(in) :: np
      INTEGER, DIMENSION(2) :: xnptobnd
      INTEGER :: ix
!
!     CHECK: DIMENSIONS WHEN SETTING FIRST BOUND
      ix = xtoix(xt)
      xnptobnd(1) = ix - INT((np-1)/2)
!
      IF (xnptobnd(1) < 0) THEN
        xnptobnd(1) = 0
      END IF
!
!     SET SECOND BOUND BY FIRST
      xnptobnd(2) = xnptobnd(1) + np -1
!
    END FUNCTION
!
!
! --------------------------------------------------------------------------- !
!
    FUNCTION ynptobnd(yt,np)
!
      INTEGER, INTENT(in) :: yt
      INTEGER, INTENT(in) :: np
      INTEGER, DIMENSION(2) :: ynptobnd
      INTEGER :: iy
!
!     CHECK: DIMENSIONS WHEN SETTING FIRST BOUND
      iy = ytoiy(yt)
      ynptobnd(1) = iy - int((np-1)/2)
!
      IF (ynptobnd(1) < 0) THEN
        ynptobnd(1) = 0
      END IF
!
!     SET SECOND BOUND BY FIRST
      ynptobnd(2) = ynptobnd(1) + np -1
!
    END FUNCTION
!
!
! --------------------------------------------------------------------------- !
!
    FUNCTION restobnd(zt,dr)
!
      INTEGER,INTENT(in) :: zt
      INTEGER, INTENT(in) :: dr
      INTEGER, DIMENSION(2) :: restobnd
      INTEGER :: iz
!
      iz = ztoiz(zt)
!
!     DEFINE SURFACE RES
      IF (dr > 0) THEN
        restobnd(1) = 0
        restobnd(2) = iz
      ELSE
        restobnd(1) = 0
        restobnd(2) = nnz
      END IF
!
    END FUNCTION
!
!
! --------------------------------------------------------------------------- !
!
    FUNCTION ztoiz(zt)
!
      INTEGER, INTENT(in) :: zt
      INTEGER :: ztoiz
!
      ztoiz = INT(zt/dz)
!
    END FUNCTION
!
!
! --------------------------------------------------------------------------- !
!
    FUNCTION xtoix(xt)
!
      INTEGER, INTENT(in) :: xt
      INTEGER :: xtoix
!
      xtoix = INT(xt/dx)
!
    END FUNCTION
!
!
! --------------------------------------------------------------------------- !
!
    FUNCTION ytoiy(yt)
!
      INTEGER, INTENT(in) :: yt
      INTEGER :: ytoiy
!
      ytoiy = int(yt/dy)
!
    END FUNCTION
!
END MODULE
