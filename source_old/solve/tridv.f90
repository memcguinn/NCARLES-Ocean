! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine describes the tridiagonal matrix solver with         !
!         multiple vectors. Nore that j and i loops are reversed from cray     !
!         version.                                                             !
!                                                                              !
!         INPUTS: n     Size of a,b,d, and r                                   !
!                 b     Below diagonal elements (b(1) not used)                !
!                 d     Diagonal elements                                      !
!                 a     Above diagonal elements (a(n) not used)                !
!                 r     Right hand side                                        !
!                 j1:j2 Range of input vectors                                 !
!                                                                              !
!         OUTPUT: r     Solution vector                                        !
! ============================================================================ !
!
SUBROUTINE tridv(b,d,a,r,n,j1,j2)
!
    REAL a(n,j1:j2), b(n,j1:j2), d(n,j1:j2), r(n,j1:j2)
!
! --------------------------------------------------------------------------- !
!
  IF (n .LE. 1) THEN
    DO j=j1,j2
      r(1,j) = r(1,j)/d(1,j)
    END DO
!
  ELSE
!
    DO j=j1,j2
      d(1,j) = 1.0/d(1,j)
    END DO
!
    DO j=j1,j2
      DO i=2,n
        fac    = b(i,j)*d(i-1,j)
        d(i,j) = 1.0/(d(i,j) - fac*a(i-1,j))
        r(i,j) = r(i,j) - fac*r(i-1,j)
      END DO
    END DO
!
    DO j=j1,j2
      r(n,j) = r(n,j)*d(n,j)
    END DO
!
    DO j=j1,j2
      DO i=n-1,1,-1
        r(i,j) = d(i,j)*(r(i,j) - a(i,j)*r(i+1,j))
      END DO
    END DO
  END IF
!
RETURN
END SUBROUTINE tridv
