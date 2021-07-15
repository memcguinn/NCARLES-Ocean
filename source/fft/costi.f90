!     SUBROUTINE COSTI(N,WSAVE)
!
!     SUBROUTINE COSTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     SUBROUTINE COST. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF SMALL PRIMES.
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             DIFFERENT WSAVE ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
!             OF N. THE CONTENTS OF WSAVE MUST NOT BE CHANGED BETWEEN
!             CALLS OF COST.
!
SUBROUTINE costi (N,WSAVE)

    DIMENSION WSAVE(*)

  PI = PIMACH(DUM)
  IF (N .LE. 3) THEN
    CONTINUE
  ELSE
    NM1 = N-1
    NP1 = N+1
    NS2 = N/2
    DT = PI/FLOAT(NM1)
    FK = 0.
    DO K=2,NS2
      KC = NP1-K
      FK = FK+1.
      WSAVE(K) = 2.*SIN(FK*DT)
      WSAVE(KC) = 2.*COS(FK*DT)
    END DO
    CALL rffti (NM1,WSAVE(N+1))
  END IF

RETURN
END