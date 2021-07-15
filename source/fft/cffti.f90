!     SUBROUTINE CFFTI(N,WSAVE)
!
!     SUBROUTINE CFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     BOTH CFFTF AND CFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4*N+15
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH CFFTF AND CFFTB
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF CFFTF OR CFFTB.
!
SUBROUTINE cffti (N,WSAVE)

    DIMENSION WSAVE(*)

  IF (N .EQ. 1) THEN
    CONTINUE
  ELSE
    IW1 = N+N+1
    IW2 = IW1+N+N
    CALL cffti1 (N,WSAVE(IW1),WSAVE(IW2))
  END IF

RETURN
END SUBROUTINE cffti