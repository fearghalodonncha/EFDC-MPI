
      SUBROUTINE VELPLTH  
C  
C CHANGE RECORD  
C  ADDED REAL FLAGS RSSBCE(L),RSSBCW(L),RSSBCN(L),RSSBCS(L)  
C  TO MODIFIED  THE OUTPUTED CELL CENTER VELOCITY FOR CELLS HAVE SOURCE/  
C **  SUBROUTINE VELPLTH WRITES A HORIZONTAL INSTANTANEOUS VELOCITY  
C **  VECTOR FILE  
C  
      USE GLOBAL 
      INTEGER*4 VER  
      INTEGER*8 TIME
       
      DIMENSION DBS(10)  
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7
      REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: UTMPS_AV
      REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: VTMPS_AV
      REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WOUT
     
      IF(.NOT.ALLOCATED(UTMPS_AV)) THEN
        ALLOCATE(UTMPS_AV(KCM))
        ALLOCATE(VTMPS_AV(KCM))
        ALLOCATE(WOUT(KCM))
      END IF
      UTMPS_AV = 0.
      VTMPS_AV = 0.
      WOUT = 0.
C  
      IF(IVPHXY.LE.2)THEN  
        IF(JSVPH.NE.1)GOTO 300  
C  
C **  WRITE HEADINGS  
C  
        TITLE1='INSTANTANEOUS HORIZ VELOCITY CM/S '  
        TITLE2='INSTANTANEOUS BOTTOM STRESS CM2/S2'  
        TITLE3='BEDLOAD TRANSPORT KG/S'  
        TITLE4='DEPTH INTEGRAED SED TRANS KG/S'  
        TITLE5='EFFECTIVE BOTTOM ROUGHNESS CM'  
        TITLE6='CURRENT SHEAR VELOCITY CM/S'  
        TITLE7='WAVE-CURRENT SHEAR VELOCITY CM/S'  
        IF(ISVPH.EQ.1) LINES1=LA-1  
        IF(ISVPH.EQ.2) LINES1=NRC  
        IF(ISVPH.EQ.3) LINES1=NBC  
        LEVELS=2  
        LEVELT=1  
        DBS(1)=0.  
        DBS(2)=99.  
        OPEN(10,FILE='VELVECH'//ans(partid2)//'.OUT',STATUS='UNKNOWN')  
        CLOSE(10,STATUS='DELETE')  
        IF (FLIU ) THEN
          OPEN(10,FILE='RUNOFF_WIND'//ans(partid2)//'.OUT',
     &             STATUS='UNKNOWN')  
          CLOSE(10,STATUS='DELETE')  
        END IF
        IF (DTFLAG) THEN  ! Output vertical velocity as well 
          OPEN(10,FILE='VELVERT'//ans(partid2)//'.OUT',STATUS='UNKNOWN')  
          CLOSE(10,STATUS='DELETE')  
        ENDIF
        IF(ISTRAN(7).GT.0)THEN  
          OPEN(10,FILE='SBLVECH.OUT',STATUS='UNKNOWN')  
          CLOSE(10,STATUS='DELETE')  
          OPEN(10,FILE='SBLVECH.OUT',STATUS='UNKNOWN')  
          WRITE (10,99) TITLE3  
          WRITE (10,101)LINES1,LEVELS  
          WRITE (10,250)(DBS(L),L=1,LEVELS)  
          CLOSE(10)  
        ENDIF  
        IF(ISWAVE.GE.1)THEN  
          OPEN(10,FILE='ZBREFFH.OUT',STATUS='UNKNOWN')  
          CLOSE(10,STATUS='DELETE')  
          OPEN(10,FILE='ZBREFFH.OUT',STATUS='UNKNOWN')  
          WRITE (10,99) TITLE5  
          WRITE (10,101)LINES1,LEVELT  
          WRITE (10,250)(DBS(L),L=1,LEVELT)  
          CLOSE(10)  
        ENDIF  
        IF(ISWAVE.GE.1)THEN  
          OPEN(10,FILE='CCUSTRH.OUT',STATUS='UNKNOWN')  
          CLOSE(10,STATUS='DELETE')  
          OPEN(10,FILE='CCUSTRH.OUT',STATUS='UNKNOWN')  
          WRITE (10,99) TITLE6  
          WRITE (10,101)LINES1,LEVELT  
          WRITE (10,250)(DBS(L),L=1,LEVELT)  
          CLOSE(10)  
        ENDIF  
        IF(ISWAVE.GE.1)THEN  
          OPEN(10,FILE='WCUSTRH.OUT',STATUS='UNKNOWN')  
          CLOSE(10,STATUS='DELETE')  
          OPEN(10,FILE='WCUSTRH.OUT',STATUS='UNKNOWN')  
          WRITE (10,99) TITLE7  
          WRITE (10,101)LINES1,LEVELT  
          WRITE (10,250)(DBS(L),L=1,LEVELT)  
          CLOSE(10)  
        ENDIF  
        JSVPH=0  
  300   CONTINUE  
        IF(ISDYNSTP.EQ.0)THEN  
          TIME=INT(DT*N)+INT(TCON*TBEGIN)  
        ELSE  
          TIME=TIMESEC
        ENDIF  
        OPEN(10,FILE='VELVECH'//ans(partid2)//'.OUT',POSITION='APPEND')  
        WRITE (10,*)N,TIME,PARTID,LA 
        IF (FLIU) THEN
          OPEN(100,FILE='RUNOFF_WIND'//ans(partid2)//'.OUT',
     &             POSITION='APPEND')  
          WRITE (100,*)N,TIME,PARTID,NQSIJ
          Write(100,*) 'Write wind stress details at this point for MH'
          WRITE(100,*) WNDVELE(LA),WNDVELN(LA),TATMT(LA) 
        END IF
        IF (DTFLAG) THEN
           OPEN(101,FILE='VELVERT'//ans(partid2)//'.OUT',
     &            POSITION='APPEND') 
           WRITE (101,*)N,TIME,PARTID,LA 
        END IF
        IF(ISTRAN(7).GT.0)THEN  
          OPEN(12,FILE='SBLVECH.OUT',POSITION='APPEND')  
          WRITE (12,100)N,TIME  
        ENDIF  
        IF(ISWAVE.GE.1)THEN  
          OPEN(14,FILE='ZBREFFH.OUT',POSITION='APPEND')  
          WRITE (14,100)N,TIME  
          OPEN(15,FILE='CCUSTRH.OUT',POSITION='APPEND')  
          WRITE (15,100)N,TIME  
          OPEN(16,FILE='WCUSTRH.OUT',POSITION='APPEND')  
          WRITE (16,100)N,TIME  
        ENDIF  
        QBOTTMP=100./CTURB3  
        IF(IVPHXY.EQ.0)THEN  
          DO L=2,LA  
            LN=LNC(L)  
            UTMPS=50.*STCUV(L)*(RSSBCE(L)*U(LEAST(L),KC)+RSSBCW(L)*U(L,KC))  
            VTMPS=50.*STCUV(L)*(RSSBCN(L)*V(LN ,KC)+RSSBCS(L)*V(L,KC))  
            VELEKC=CUE(L)*UTMPS+CVE(L)*VTMPS  
            VELNKC=CUN(L)*UTMPS+CVN(L)*VTMPS  
            UTMPB=50.*STCUV(L)*(RSSBCE(L)*U(LEAST(L),1)+RSSBCW(L)*U(L,1))  
            VTMPB=50.*STCUV(L)*(RSSBCN(L)*V(LN ,1)+RSSBCS(L)*V(L,1))  
            VELEKB=CUE(L)*UTMPB+CVE(L)*VTMPB  
            VELNKB=CUN(L)*UTMPB+CVN(L)*VTMPB  
            UTMPA=50.*STCUV(L)*(RSSBCE(L)*UHE(LEAST(L))*HUI(LEAST(L))  
     &          +RSSBCW(L)*UHE(L)*HUI(L))  
            VTMPA=50.*STCUV(L)*(RSSBCN(L)*VHE(LN )*HVI(LN )  
     &          +RSSBCS(L)*VHE(L)*HVI(L))  
            TUTMP=5000.*STCUV(L)*(RSSBCE(L)*TBX(LEAST(L))+RSSBCW(L)*TBX(L))  
            TVTMP=5000.*STCUV(L)*(RSSBCN(L)*TBY(LN )+RSSBCS(L)*TBY(L))  
            VELEAV=CUE(L)*UTMPA+CVE(L)*VTMPA  
            VELNAV=CUN(L)*UTMPA+CVN(L)*VTMPA  
C  
C WRITE VELVECH.OUT  
C  
            IF(KC.GT.1)WRITE(10,201)  
     &          VELEKC,VELNKC,VELEKB,VELNKB,VELEAV,VELNAV  
            IF(KC.EQ.1)WRITE(10,200)IL(L),JL(L),VELEKB,VELNKB  
C  
C WRITE VELVECH.COC  
C  
            UTMP=5000.*STCUV(L)*(RSSBCE(L)*TBX(LEAST(L))+RSSBCW(L)*TBX(L))  
            VTMP=5000.*STCUV(L)*(RSSBCN(L)*TBY(LN )+RSSBCS(L)*TBY(L))  
            VELEKC=CUE(L)*UTMP+CVE(L)*VTMP  
            VELNKC=CUN(L)*UTMP+CVN(L)*VTMP  
            TMPV1=10000.*TAUB(L)  
            TMPV2=10000.*TAUBSED(L)  
            TMPV3=10000.*TAUBSND(L)  
C  
C WRITE TAUVECH.OUT  
C  
            IF(ISTRAN(7).GE.1) THEN  
              UTMP=0.0005*STCUV(L)*(RSSBCE(L)*QSBDLDX(LEAST(L),1)  
     &            +RSSBCW(L)*QSBDLDX(L  ,1))  
              VTMP=0.0005*STCUV(L)*(RSSBCN(L)*QSBDLDY(LN ,1)  
     &            +RSSBCS(L)*QSBDLDY(L  ,1))  
              VELEKC=CUE(L)*UTMP+CVE(L)*VTMP  
              VELNKC=CUN(L)*UTMP+CVN(L)*VTMP  
              UTMP=0.0005*STCUV(L)*(RSSBCE(L)*QSBDLDX(LEAST(L),2)  
     &            +RSSBCW(L)*QSBDLDX(L,2))  
              VTMP=0.0005*STCUV(L)*(RSSBCN(L)*QSBDLDY(LN ,2)  
     &            +RSSBCS(L)*QSBDLDY(L,2))  
              VELEKB=CUE(L)*UTMP+CVE(L)*VTMP  
              VELNKB=CUN(L)*UTMP+CVN(L)*VTMP  
              WRITE(12,200)VELEKC,VELNKC,  
     &            VELEKB,VELNKB  
            END IF  
            IF(ISWAVE.EQ.2)THEN  
              ZBREFF=100.*ZBRE(L)  
              WRITE(14,201)ZBREFF  
              QTURBC=QBOTTMP*QQSQR(L,0)  
              WRITE(15,201)QTURBC  
              QTURBC=QBOTTMP*QQWV2(L)  
              WRITE(16,201)QTURBC  
            ENDIF  
          ENDDO  
        ENDIF  
        IF(IVPHXY.EQ.1)THEN  
          DO L=2,LA  
            LE=LEAST(L)
            LN=LNC(L)  
            UTMPS=50.*STCUV(L)*(RSSBCE(L)*U(LEAST(L),KC)+RSSBCW(L)*U(L,KC))  
            VTMPS=50.*STCUV(L)*(RSSBCN(L)*V(LN ,KC)+RSSBCS(L)*V(L,KC))  
            VELEKC=CUE(L)*UTMPS+CVE(L)*VTMPS  
            VELNKC=CUN(L)*UTMPS+CVN(L)*VTMPS  
            UTMPB=50.*STCUV(L)*(RSSBCE(L)*U(LEAST(L),1)+RSSBCW(L)*U(L,1))  
            VTMPB=50.*STCUV(L)*(RSSBCN(L)*V(LN ,1)+RSSBCS(L)*V(L,1))  
            VELEKB=CUE(L)*UTMPB+CVE(L)*VTMPB  
            VELNKB=CUN(L)*UTMPB+CVN(L)*VTMPB  
            UTMPA=50.*STCUV(L)*(RSSBCE(L)*UHE(LEAST(L))*HUI(LEAST(L))  
     &          +RSSBCW(L)*UHE(L)*HUI(L))  
            VTMPA=50.*STCUV(L)*(RSSBCN(L)*VHE(LN )*HVI(LN )  
     &          +RSSBCS(L)*VHE(L)*HVI(L))  
            TUTMP=5000.*STCUV(L)*(RSSBCE(L)*TBX(LEAST(L))+RSSBCW(L)*TBX(L))  
            TVTMP=5000.*STCUV(L)*(RSSBCN(L)*TBY(LN )+RSSBCS(L)*TBY(L))  
            VELEAV=CUE(L)*UTMPA+CVE(L)*VTMPA  
            VELNAV=CUN(L)*UTMPA+CVN(L)*VTMPA  
C  
C WRITE VELVECH.OUT  
C  
            do k = 1,kc
              UTMPS_AV(K) = 0.5*STCUV(L)*(RSSBCE(L)*U(LE,K)+
     &                      RSSBCW(L)*U(L,K) )
              VTMPS_AV(K) = 0.5*STCUV(L)*(RSSBCN(L)*V(LN,K)+
     &                      RSSBCS(L)*V(L,K) )
              WOUT(K)  = 0.5*(W(L,K) + W(L,K-1)) ! Average as per Arakawa grid

            end do
            I = IL(L)
            J = JL(L)
            IF(KC.GT.1)WRITE(10,210)XPAR(I),YPAR(J),
     &        (utmps_av(K),K=1,KC),(vtmps_av(K),K=1,KC),(WOUT(K),K=1,KC)
            IF(KC.EQ.1)WRITE(10,200)XPAR(I),YPAR(J),VELEKB,VELNKB,W(L,1)
            
C  
C WRITE VELVECH.COC  
C  
            UTMP=5000.*STCUV(L)*(RSSBCE(L)*TBX(LEAST(L))+RSSBCW(L)*TBX(L))  
            VTMP=5000.*STCUV(L)*(RSSBCN(L)*TBY(LN )+RSSBCS(L)*TBY(L))  
            VELEKC=CUE(L)*UTMP+CVE(L)*VTMP  
            VELNKC=CUN(L)*UTMP+CVN(L)*VTMP  
            TMPV1=10000.*TAUB(L)  
            TMPV2=10000.*TAUBSED(L)  
            TMPV3=10000.*TAUBSND(L)  
C  
            IF(ISTRAN(7).GE.1) THEN  
              UTMP=0.0005*STCUV(L)*(RSSBCE(L)*QSBDLDX(LEAST(L),1)  
     &            +RSSBCW(L)*QSBDLDX(L  ,1))  
              VTMP=0.0005*STCUV(L)*(RSSBCN(L)*QSBDLDY(LN ,1)  
     &            +RSSBCS(L)*QSBDLDY(L  ,1))  
              VELEKC=CUE(L)*UTMP+CVE(L)*VTMP  
              VELNKC=CUN(L)*UTMP+CVN(L)*VTMP  
              UTMP=0.0005*STCUV(L)*(RSSBCE(L)*QSBDLDX(LEAST(L),2)  
     &            +RSSBCW(L)*QSBDLDX(L,2))  
              VTMP=0.0005*STCUV(L)*(RSSBCN(L)*QSBDLDY(LN ,2)  
     &            +RSSBCS(L)*QSBDLDY(L,2))  
              VELEKB=CUE(L)*UTMP+CVE(L)*VTMP  
              VELNKB=CUN(L)*UTMP+CVN(L)*VTMP  
              WRITE(12,200)IL(L),JL(L),VELEKC,VELNKC,  
     &            VELEKB,VELNKB  
            END IF  
            IF(ISWAVE.EQ.2)THEN  
              ZBREFF=100.*ZBRE(L)  
              WRITE(14,200)IL(L),JL(L),DLON(L),DLAT(L),ZBREFF  
              QTURBC=QBOTTMP*QQSQR(L,0)  
              WRITE(15,200)IL(L),JL(L),DLON(L),DLAT(L),QTURBC  
              QTURBC=QBOTTMP*QQWV2(L)  
              WRITE(16,200)IL(L),JL(L),DLON(L),DLAT(L),QTURBC  
            ENDIF  
          ENDDO  
        ENDIF 
       !IF(IVPHXY.EQ.2)THEN  
        !END IF  
        CLOSE(10)  
        IF(ISTRAN(7).GT.0)CLOSE(12)  
        CLOSE(13)  
        CLOSE(14)  
        CLOSE(15)  
        CLOSE(16)  
        CLOSE(20)  
      ENDIF  
C  
C *** EE BEGIN BLOCK  
C *** OUTPUT EFDC EXPLORER FORMAT.  DO NOT CHANGE OUTPUTS!  
C ***                               MUST EXACTLY MATCH EFDC_EXPLORER INP  
C  
      IF(IVPHXY.EQ.3)THEN  
        IF(JSVPH.EQ.1)THEN  
          LINES=LA-1  
          OPEN(10,FILE='EE_VEL.OUT',STATUS='UNKNOWN',  
     &        ACCESS='SEQUENTIAL',FORM='UNFORMATTED')  
          VER=103  
          WRITE(10)VER,IC,JC,KC,LINES  

          WRITE(10)RSSBCE,RSSBCW,RSSBCS,RSSBCN

          CLOSE(10)  
          JSVPH=0  
        ENDIF  
        IF(ISDYNSTP.EQ.0)THEN  
          TIME=INT(DT*FLOAT(N)+TCON*TBEGIN)  
        ELSE  
          TIME=TIMESEC  
        ENDIF  
        TIME=TIME/86400.  
        IF(ISDYNSTP.EQ.0)THEN  
          DELT=DT  
        ELSE  
          DELT=DTDYN  
        ENDIF  
        OPEN(10,FILE='EE_VEL.OUT',POSITION='APPEND',STATUS='OLD',  
     &      FORM='UNFORMATTED')  
        WRITE (10)N,TIME,DELT  

        ! *** Write the UVW Instantaneous Velocity Field (Unrotated)
        DO L=2,LA  
          WRITE(10)(U(L,K),V(L,K),W(L,K),K=1,KC)  
        ENDDO  

CFERG        CALL FLUSH(10)
        CLOSE(10)  
      ENDIF  

! Write river data to file for Mike H
      IF (FLIU) THEN
      DO LL=1,NQSIJ  
        NS=NQSERQ(LL)  
        L=LQS(LL)  
        I = IL(L)
        J = JL(L)
        QOUT = SUM(QSERT(:,NS))
            write(100,120) XPAR(I),YPAR(J),QOUT
       END DO
       END IF
       CLOSE(100)
       CLOSE(101)

C  
C *** EE END BLOCK  
C  
   99 FORMAT(A80)  
  100 FORMAT(I10,F12.4)  
  120 FORMAT(2I5,F12.4)  

 111  FORMAT(I10,F12.4,I5,I10)  
  101 FORMAT(2I10)  
  200 FORMAT(2I5,1X,20E14.6)  
  201 FORMAT(8E14.6)  
  210 FORMAT(2I5,1X,200E14.6)
  250 FORMAT(12E12.4)  
      RETURN  
      END  

