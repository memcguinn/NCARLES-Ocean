! ============================================================================ !
! ABOUT:                                                                       !
!         This subroutine consists of the IBM range finder to balance load.    !
! ============================================================================ !
!
SUBROUTINE range(n1,n2,nprocs,irank,ista,iend)
!
! --------------------------------------------------------------------------- !
!
  iwork1 = (n2 - n1 + 1)/nprocs
  iwork2 = MOD(n2 - n1 +1, nprocs)
  ista = irank*iwork1 + n1 + MIN(irank,iwork2)
  iend = ista + iwork1 - 1
!
  IF (iwork2 .GT. irank) iend = iend + 1
!
RETURN
!
END SUBROUTINE range
