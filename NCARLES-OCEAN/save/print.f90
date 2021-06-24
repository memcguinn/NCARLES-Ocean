! ============================================================================ !
! ABOUT:                                                                       !
! ============================================================================ !
!
SUBROUTINE print(lu,it,iz_strt,iz_end)
!
    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
!
! --------------------------------------------------------------------------- !
!
  WRITE (lu,4000)
  WRITE (lu,4100) it,time,dt,zi,tsfcc(1),uwsfc,vwsfc,wtsfc(1),   &
                  zol,hol,ucfl, vcfl, wcfl,t_ref, u_10, v_10,    &
                  cd_10,tau_x, tau_y, cpou10,stokess,   &
                  fcor, fcor_h
  WRITE (lu,4200)
!
  DO iz=iz_end,iz_strt,-1
    WRITE (lu,4300) iz,txym(iz,1)-t_ref,divz(iz),englez(iz),     &
                    eavg(iz),wtle(iz,1),wtsb(iz,1),shrz(iz),buyz(iz)
  END DO
!
  WRITE (lu,4400) tsfcc(1),wtsfc(1)
  WRITE (lu,4500) stokess,udrift,vdrift
  WRITE (lu,4600) (iz,uxym(iz)+ugal,vxym(iz),uwle(iz),uwsb(iz),  &
                  vwle(iz),vwsb(iz),iz=iz_strt,iz_end)
  WRITE (lu,4800) xksurf, nmatch, viscon, vise
  WRITE (lu,4700) (iz,dfac(iz),iz=iz_strt,iz_end)
  WRITE (lu,5100) (iz,txym(iz,2),txym(iz,3),txym(iz,4),          &
                   txym(iz,5),txym(iz,6),txym(iz,7),iz=iz_strt,iz_end)
!
RETURN
!
! ---------------------------------- FORMAT ---------------------------------- !
4000 FORMAT (30X,' --- SOLUTION ---')
4100 FORMAT (' IT=',I7,5x,'TIME (s) = ',e15.8,',  DT(s) = ',e15.6,/,    &
             10x,'ZTOP = ',e15.6,',  TSFC = ',e15.6,',  UW = ',         &
             e15.6,',  VW = ',e15.6,/,10x,'WT = ',e15.6,',  ZL =',      &
             e15.6,',  HL = ',e15.6,/,10x,'U_cfl = ',e15.6,             &
             ',  V_cfl = ',e15.6,',  W_cfl = ',e15.6,/,10x,'Theta Ref = ', &
             e15.6,/,10x,'U_10 = ',e15.6,',  V_10 = ',e15.6,            &
             ',  Cd_10 = ',e15.6,/,10x,'Tau_x = ',e15.6,',  Tau_y = ',  &
             e15.6,/,10x,'La_t = ',e15.6,',  Wave age = ',e15.6,/,10x,  &
             'Stokes amp = ',e15.6,',  Stokes wave number = ',e15.6,/,  &
             10x,'f cor = ',e15.6,',  f cor_h = ',e15.6)
!
4200 FORMAT (//,20x,'--------- HORIZONTAL MEAN VALUES ---------- ',     &
             //,2x,'IZ',4x,'T_MEAN',7x,'DIVG',8X,'LE_KE',6X,'SGS_KE',   &
             7X,'LE_WT',6X,'SGS_WT',7X,'SHRZ',8X,'BUOY')
4300 FORMAT (1X,I3,e12.4,7e12.4)
4400 FORMAT ('  SURFACE VALUE: TXYM=',F8.2,'               WTSB=',E9.2)
4500 FORMAT (/,' STOKESS = ',e12.4,' UDRIFT = ',e12.4,' VDRIFT = ',e12.4)
4600 FORMAT (//,' IZ',5x,' UXYM + UGAL',8x,' VXYM',10x,' UWLE',10x,     &
             ' UWSB',10x,' VWLE',10x,' VWSB',/,(1x,i4,6(3x,e13.6)))
4800 FORMAT (//,' XKSURF = ',e15.6,' NMATCH = ',i4,/,' VISCON = ',      &
             e15.6,' VISE = ',e15.6)
4700 FORMAT (//,'   IZ',5x,'  DFAC',/,(1x,i4,3x,e15.6))
5100 FORMAT (//,' IZ',5x,' SCALAR-1 MEAN',8x,'2-MEAN',10x,'3-MEAN',10x, &
             '4-MEAN', 10x, '5-MEAN',10x,'6-MEAN',10x,/,(1x,i4,6(3x,e13.6)))
!------------------------------------------------------------------------------!
!
END SUBROUTINE print
