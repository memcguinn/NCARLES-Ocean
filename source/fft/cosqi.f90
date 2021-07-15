!     SUBROUTINE COSQI(N,WSAVE)
!
!     SUBROUTINE COSQI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     BOTH COSQF AND COSQB. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE ARRAY TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH COSQF AND COSQB
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF COSQF OR COSQB.
!
SUBROUTINE cosqi (N,WSAVE)

    DIMENSION WSAVE(*)

  PIH = 0.5*PIMACH(DUM)
  DT = PIH/FLOAT(N)
  FK = 0.
  DO K=1,N
    FK = FK+1.
    WSAVE(K) = COS(FK*DT)
  END DO
  CALL rffti (N,WSAVE(N+1))

RETURN
END SUBROUTINE cosqi