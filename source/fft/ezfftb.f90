!     SUBROUTINE EZFFTB(N,R,AZERO,A,B,WSAVE)
!
!     SUBROUTINE EZFFTB COMPUTES A REAL PERODIC SEQUENCE FROM ITS
!     FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS
!     DEFINED BELOW AT OUTPUT PARAMETER R. EZFFTB IS A SIMPLIFIED
!     BUT SLOWER VERSION OF RFFTB.
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE OUTPUT ARRAY R.  THE METHOD IS MOST
!             EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
!
!     AZERO   THE CONSTANT FOURIER COEFFICIENT
!
!     A,B     ARRAYS WHICH CONTAIN THE REMAINING FOURIER COEFFICIENTS
!             THESE ARRAYS ARE NOT DESTROYED.
!
!             THE LENGTH OF THESE ARRAYS DEPENDS ON WHETHER N IS EVEN OR
!             ODD.
!
!             IF N IS EVEN N/2    LOCATIONS ARE REQUIRED
!             IF N IS ODD (N-1)/2 LOCATIONS ARE REQUIRED
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             IN THE PROGRAM THAT CALLS EZFFTB. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE EZFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY EZFFTF AND EZFFTB.
!
!
!     OUTPUT PARAMETERS
!
!     R       IF N IS EVEN DEFINE KMAX=N/2
!             IF N IS ODD  DEFINE KMAX=(N-1)/2
!
!             THEN FOR I=1,...,N
!
!                  R(I)=AZERO PLUS THE SUM FROM K=1 TO K=KMAX OF
!
!                  A(K)*COS(K*(I-1)*2*PI/N)+B(K)*SIN(K*(I-1)*2*PI/N)
!
!     ********************* COMPLEX NOTATION **************************
!
!             FOR J=1,...,N
!
!             R(J) EQUALS THE SUM FROM K=-KMAX TO K=KMAX OF
!
!                  C(K)*EXP(I*K*(J-1)*2*PI/N)
!
!             WHERE
!
!                  C(K) = .5*CMPLX(A(K),-B(K))   FOR K=1,...,KMAX
!
!                  C(-K) = CONJG(C(K))
!
!                  C(0) = AZERO
!
!                       AND I=SQRT(-1)
!
!     *************** AMPLITUDE - PHASE NOTATION ***********************
!
!             FOR I=1,...,N
!
!             R(I) EQUALS AZERO PLUS THE SUM FROM K=1 TO K=KMAX OF
!
!                  ALPHA(K)*COS(K*(I-1)*2*PI/N+BETA(K))
!
!             WHERE
!
!                  ALPHA(K) = SQRT(A(K)*A(K)+B(K)*B(K))
!
!                  COS(BETA(K))=A(K)/ALPHA(K)
!
!                  SIN(BETA(K))=-B(K)/ALPHA(K)
!
SUBROUTINE ezfftb (N,R,AZERO,A,B,WSAVE)

    DIMENSION R(*), A(*), B(*), WSAVE(*)

  IF (N-2 .LT. 0) THEN
    R(1) = AZERO
  ELSE IF (N-2 .EQ. 0) THEN
    R(1) = AZERO+A(1)
    R(2) = AZERO-A(1)
  ELSE !(N-2 .GT. 0)
    NS2 = (N-1)/2
    DO I=1,NS2
      R(2*I) = .5*A(I)
      R(2*I+1) = -.5*B(I)
    END DO
    R(1) = AZERO
    IF (MOD(N,2) .EQ. 0) R(N) = A(NS2+1)
    CALL rfftb (N,R,WSAVE(N+1))
  END IF

RETURN
END SUBROUTINE ezfftb
