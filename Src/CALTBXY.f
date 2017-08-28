      SUBROUTINE CALTBXY(ISTL_,IS2TL_)  
!  
! **  SUBROUTINE CALTBXY CALCULATES BOTTOM FRICTION OR DRAG  
! **  COEFFICIENTS IN QUADRATIC LAW FORM REFERENCED TO NEAR  
! **  BOTTOM OR DEPTH AVERAGED HORIZONTAL VELOCITIES  
! **  FOR VEGETATION RESISTANCE IN DEPTH INTEGRATED FLOW  
! **  THE COEFFICIENT REPRESENTS BOTTOM AND WATER COLUMN VEGETATION  
! **  RESISTANCE  
! CHANGE RECORD  
!  REMOVED DRAG COEFFICIENT CONSTRAINT FOR MULIPLE LAYER ROUGH  
!   BOUNDARIES WHEN DYNAMIC TIME STEPPING IS ACTIVE  
!  FIXED POSSIBLE DIVIDE BY ZERO FOR SUB GRID CHANNEL FRICTION IN  
!  ABSENCE OF VEGETATION RESISTANCE  
!  ADDED DRY CELL BYPASS AND CONSISTENT INITIALIZATION OF DRY VALUES  
!  
      USE GLOBAL  
      USE OMP_LIB
	IMPLICIT NONE
      DOUBLE PRECISION LST,LEND,foo
	INTEGER::ISTL_,IS2TL_,L,K,LS,M,LW,LE,LN,LNW,LSE,MW,MS
	INTEGER::NMD,LHOST,LCHNU,LCHNV,MH,MU,MV,NTMP
	INTEGER::LZBMIN,LCDMAX,LCDMIN,LZBMAX,JWCBLV,JWCBLU
	REAL::CDLIMIT,CDTOTUM,CDTOTVM,CDMAXUM,CDMAXVM
	REAL::ZBRATU,ZBRATV,UMAGTMP,VMAGTMP,CDMAXU,CDMAXV
	REAL::HURTMP,HVRTMP,HUDZBR,HVDZBR,VTMPATU,UTMPATV,CPVEGU,RVEGUM
	REAL::CPVEGV,RVEGVM,HVGTC,HVGTW,HVGTS,VISEXP,VISFAC,VISMUDU
	REAL::VISMUDV,SEDTMP,CSEDVIS,VISDHU,VISDHV,DZHUDZBR,DZHVDZBR
	REAL::FRACLAY,FHLAYC,FHLAYW,FHLAYS,WCHAN,RLCHN,HCHAN,STBXCH
	REAL::FXVEGCH,STBYCH,FYVEGCH,TMPVALW,WVFACT,QQWCTMP,TWCTMP
	REAL::AEXTMP,TMPVAL,USTARC,CDRGTMP,TAUBTMP,TAUE,RIPAMP
	REAL::RIPSTP,RIPFAC,ZBRMAX,ZBRMIN,CDRGMAX,ZBREU
	REAL::CDRGMIN,WVDTMP,RKZTURB,UTMP,VTMP,DWVDZ,DWUDZ,DWVD2Z
	REAL::DWUD2Z,HZRVDZ,HZRUDZ,ZDHZRV,ZDHZRU,ZBREV,HZREFV,HZREFU
	REAL::QWDQCV,QWDQCU,QCTMPV,QCTMPU,HOTLYMN,HOTLYMX,CDTMPVY
	REAL::BOTTMP,DWVDHR,DWUDHR,QWCTMPV,QWCTMPU
	REAL::CDTMPV,CDTMPU,COSWC,CURANG,CDTMPUX
	REAL::WVDELV,WVDELU,TAUTMP

      DELT=DT2  
      ISUD=1  
      IF(ISTL_.NE.3)THEN  
        DELT=DT  
        ISUD=0  
      ENDIF  
      IF(IS2TL_.EQ.1)THEN  
        IF(ISDYNSTP.EQ.0)THEN  
          DELT=DT  
        ELSE  
          DELT=DTDYN  
        END IF  
        ISUD=1  
      ENDIF  
      DELTI=1./DELT  
!  
! **  IF WAVE-CURRENT BBL MODEL IS ACTIVE, GOTO WAVE CURRENT BBL  
!  
      IF(ISWCBL.GE.1) GOTO 1947  
!  
! **  INITIALIZE IMPLICIT BOTTOM FRICTION AND SET DIAGNOSTIC FILES  
! **  ON FIRST CALL  
!  
      IF(JSTBXY.EQ.1) GOTO 100  
      IF(ISITB.GE.1)THEN  
        IF(ISITB.EQ.1)THEN  
          RITB1=0.45  
          RITB=0.55  
          CDLIMIT=1.  
        ELSE  
          RITB1=0.0  
          RITB=1.0  
          CDLIMIT=10.  
        ENDIF  
      ELSE  
        RITB1=1.0  
        RITB=0.0  
        CDLIMIT=0.5  
      ENDIF  
      IF(ISVEG.GE.2)THEN  
        OPEN(1,FILE='CBOT.LOG',STATUS='UNKNOWN')  
        CLOSE(1,STATUS='DELETE')  
      ENDIF  
      DO L=2,LA  
        STBXO(L)=STBX(L)  
        STBYO(L)=STBY(L)  
      ENDDO  
      DO L=1,LC  
        STBX(L)=0.  
        STBY(L)=0.  
      ENDDO  
      DO K=1,KC  
        DO L=1,LC  
          FXVEG(L,K)=0.  
          FYVEG(L,K)=0.  
        ENDDO  
      ENDDO  
      N=-2  
      JSTBXY=1  
  100 CONTINUE  
      IF(ISITB.GE.1)THEN  
        IF(ISITB.EQ.1)THEN  
          CDLIMIT=10.  
        ELSE  
          CDLIMIT=100.  
        ENDIF  
      ELSE  
        CDLIMIT=0.5  
      ENDIF  
!  
! **  INITIALIZED DIAGNOSTICS FOR STANDARD AND VEGE  
! **  RESISTANCE CALCULATION  
!  
      IF(ISVEG.GE.2)THEN  
        OPEN(1,FILE='CBOT.LOG',POSITION='APPEND',STATUS='UNKNOWN')  
      ENDIF  
      CDTOTUM=0.  
      CDTOTVM=0.  
      CDMAXUM=0.  
      CDMAXVM=0.  
      IF(ISVEG.EQ.0) UVEGSCL=1.E-12  
      IF(KC.GT.1) GOTO 200  
!  
! **  NORMAL ENTRY INTO STANDARD AND VEGE RESISTANCE CALCULATION  
! **  FOR SINGLE LAYER  
! **  VEGETATION DRAG  
!           CALCULATE R FOR LAMINAR FLOW  
!           CALCULATE R FOR LAMINAR FLOW  
! **  END VEGETATION DRAG  
! **  NORMAL ENTRY INTO STANDARD AND VEGE RESISTANCE CALCULATION  
! **  FOR SINGLE LAYER  
!  
!$OMP PARALLEL PRIVATE(LS,LW,ZBRATU,ZBRATV,
!$OMP& UMAGTMP,VMAGTMP,CDMAXU,CDMAXV,HURTMP,HVRTMP,
!$OMP& HUDZBR,HVDZBR)
!$OMP DO SCHEDULE(static,chunksize)
      DO L=2,LA  
        IF(LMASKDRY(L))THEN  
          LS=LSC(L)  
          LW=LWEST(L)
          ZBRATU=0.5*(DXP(LW)*ZBR(LW)+DXP(L)*ZBR(L))*DXIU(L)  
          ZBRATV=0.5*(DYP(LS )*ZBR(LS )+DYP(L)*ZBR(L))*DYIV(L)  
          UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L)+1.E-12 )  
          VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1)+1.E-12 )  
          CDMAXU=CDLIMIT*STBXO(L)*H1U(L)/( DELT*UMAGTMP )  
          CDMAXV=CDLIMIT*STBYO(L)*H1V(L)/( DELT*VMAGTMP )  
          HURTMP=MAX(ZBRATU,H1U(L))  
          HVRTMP=MAX(ZBRATV,H1V(L))  
          HUDZBR=HURTMP/ZBRATU  
          IF(HUDZBR.LT.7.5) HUDZBR=7.5  
          HVDZBR=HVRTMP/ZBRATV  
          IF(HVDZBR.LT.7.5) HVDZBR=7.5  
          STBX(L)=STBXO(L)*.16/( (LOG( HUDZBR ) -1.)**2)  
          STBY(L)=STBYO(L)*.16/( (LOG( HVDZBR ) -1.)**2)  
          STBX(L)=MIN(CDMAXU,STBX(L))  
          STBY(L)=MIN(CDMAXV,STBY(L))  
        ENDIF  
      ENDDO  
!$OMP END PARALLEL
      IF(ISVEG.GE.1)THEN  
        K=1  
        DO L=2,LA  
          IF(LMASKDRY(L))THEN  
            M=MVEGL(L)  
            FXVEG(L,K)=0.  
            FYVEG(L,K)=0.  
!  
! *** DSLLC BEGIN BLOCK  
!  
            IF(M.NE.MVEGOW.AND.M.NE.0.AND.M.LT.91)THEN  
              LW=LWEST(L)  
              LE=LEAST(L)  
              LS=LSC(L)  
              LN=LNC(L)  
              LNW=LNWC(L)  
              LSE=LSEC(L)  
              MW=MVEGL(LW)  
              MS=MVEGL(LS)  
              VTMPATU=0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))  
              UTMPATV=0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))  
              UMAGTMP=SQRT( U(L,K)*U(L,K)+VTMPATU*VTMPATU +1.E-12 )  
              VMAGTMP=SQRT( UTMPATV*UTMPATV+V(L,K)*V(L,K) +1.E-12 )  
              UMAGTMP=MAX(UMAGTMP,UVEGSCL)  
              VMAGTMP=MAX(VMAGTMP,UVEGSCL)  
              CDMAXU=CDLIMIT*STBXO(L)*H1U(L)/( DELT*UMAGTMP )  
              CDMAXV=CDLIMIT*STBYO(L)*H1V(L)/( DELT*VMAGTMP )  
              IF(N.EQ.-2)THEN  
                VTMPATU=0.25*(V1(L,K)+V1(LW,K)+V1(LN,K)+V1(LNW,K))  
                UTMPATV=0.25*(U1(L,K)+U1(LE,K)+U1(LS,K)+U1(LSE,K))  
                UMAGTMP=SQRT( U1(L,K)*U1(L,K)+VTMPATU*VTMPATU+1.E-12 )  
                VMAGTMP=SQRT( UTMPATV*UTMPATV+V1(L,K)*V1(L,K)+1.E-12 )  
              ENDIF  
CJH           CPVEGU=0.5 ! CHANGED DEFINITION
              CPVEGU=1.0
              IF(ISVEGL.EQ.1) CPVEGU=CPVEGU + 10.E-6/(  
     &            (BPVEG(MW)+BPVEG(M))*UMAGTMP )  
              IF(CPVEGU.GT.1.0)THEN  
!  
!            CALCULATE R FOR LAMINAR FLOW  
!  
                CPVEGU=CPVEGU-0.5  
                RVEGUM=0.  
              ENDIF  
              CPVEGU=SCVEG(M)*CPVEGU  
CHH           CPVEGV=0.5 ! CHANGED DEFINITION
              CPVEGV=1.0
              IF(ISVEGL.EQ.1) CPVEGV=CPVEGV + 10.E-6/(  
     &            (BPVEG(MS)+BPVEG(M))*VMAGTMP )  
              IF(CPVEGV.GT.1.0)THEN  
!  
!            CALCULATE R FOR LAMINAR FLOW  
!  
                CPVEGV=CPVEGV-0.5  
                RVEGVM=0.  
              ENDIF  
              CPVEGV=SCVEG(M)*CPVEGV  
              HVGTC=MIN(HPVEG(M),HP(L))  
              HVGTW=MIN(HPVEG(MW),HP(LWEST(L)))  
              HVGTS=MIN(HPVEG(MS),HP(LS))  
              FXVEG(L,K)=0.25*CPVEGU*( DXP(L)*(BDLPSQ(M)*HVGTC/PVEGZ(M))  
     &            +DXP(LWEST(L))*(BDLPSQ(MW)*HVGTW/PVEGZ(MW)) )*DXIU(L)  
              FYVEG(L,K)=0.25*CPVEGV*( DYP(L)*(BDLPSQ(M)*HVGTC/PVEGZ(M))  
     &            +DYP(LS)*(BDLPSQ(MS)*HVGTS/PVEGZ(MS)) )*DYIV(L)  
              FXVEG(L,K)=MIN(FXVEG(L,K),CDMAXU)  
              FYVEG(L,K)=MIN(FYVEG(L,K),CDMAXU)
          ENDIF  
!  
! *** DSLLC END BLOCK  
!  
          ENDIF  
        ENDDO  
      ENDIF  
      GOTO 300  
!  
! **  NORMAL ENTRY INTO STANDARD AND VEGE RESISTANCE CALCULATION  
! **  FOR MULTIPLE LAYER  
!  
  200 CONTINUE  
!  
! **  BEGIN SMOOTH DRAG FORMULATION  
!  
      VISEXP=2./7.
      VISFAC=0.0258*(COEFTSBL**VISEXP)
!
      DO L=2,LA  
        IF(LMASKDRY(L))THEN  
          IF(ZBR(L).LE.1.E-6)THEN  
            UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L)+1.E-12 )  
            VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1)+1.E-12 )  
            CDMAXU=CDLIMIT*STBXO(L)*H1U(L)/( DELT*UMAGTMP )  
            CDMAXV=CDLIMIT*STBYO(L)*H1V(L)/( DELT*VMAGTMP )  
            VISMUDU=VISMUD
            VISMUDV=VISMUD
            IF(ISMUD.GE.1)THEN  
              SEDTMP=0.5*(SED(L,1,1)+SED(LWEST(L),1,1))  
              VISMUDU=CSEDVIS(SEDTMP)  
              SEDTMP=0.5*(SED(L,1,1)+SED(LSC(L),1,1))  
              VISMUDV=CSEDVIS(SEDTMP)  
            ENDIF  
! **  DELETED COMMENTED OUT LINES & UNUSED VARIABLES
            VISDHU=0.0
            VISDHV=0.0
            IF(UMAGTMP.GT.0.0) VISDHU=(VISMUDU*HUI(L)/UMAGTMP)*VISEXP
            IF(VMAGTMP.GT.0.0) VISDHV=(VISMUDV*HVI(L)/VMAGTMP)*VISEXP
            STBX(L)=VISFAC*AVCON*STBXO(L)*VISDHU
            STBY(L)=VISFAC*AVCON*STBYO(L)*VISDHV
            STBX(L)=MIN(CDMAXU,STBX(L))  
            STBY(L)=MIN(CDMAXV,STBY(L))  
          ENDIF  
        ENDIF  
      ENDDO  
!  
! **  END SMOOTH DRAG FORMULATION  
!
! **  BEGIN ROUGH DRAG FORMULATION  
!  
      IF(N.NE.-2)THEN
!$OMP PARALLEL PRIVATE(LS,LW,ZBRATU,ZBRATV,
!$OMP& UMAGTMP,VMAGTMP,CDMAXU,CDMAXV,HURTMP,HVRTMP,
!$OMP& DZHUDZBR,DZHVDZBR)
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)

	    DO L=2,LA  
          IF(LMASKDRY(L))THEN  
            LS=LSC(L)  
            LW=LWEST(L)
            IF(ZBR(L).GT.1.E-6)THEN  
              ZBRATU=0.5*(DXP(LW)*ZBR(LW)+DXP(L)*ZBR(L))*DXIU(L)  
              ZBRATV=0.5*(DYP(LS )*ZBR(LS )+DYP(L)*ZBR(L))*DYIV(L)  
              UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L)+1.E-12 )  
              VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1)+1.E-12 )  
              CDMAXU=CDLIMIT*STBXO(L)*H1U(L)/( DELT*UMAGTMP )  
              CDMAXV=CDLIMIT*STBYO(L)*H1V(L)/( DELT*VMAGTMP )  
              HURTMP=MAX(ZBRATU,H1U(L))  
              HVRTMP=MAX(ZBRATV,H1V(L))  
              DZHUDZBR=1.+0.5*DZC(1)*HURTMP/ZBRATU  
              DZHVDZBR=1.+0.5*DZC(1)*HVRTMP/ZBRATV  
!  
              STBX(L)=AVCON*STBXO(L)*.16/((LOG(DZHUDZBR))**2)  
              STBY(L)=AVCON*STBYO(L)*.16/((LOG(DZHVDZBR))**2)  
              STBX(L)=MIN(CDMAXU,STBX(L))  
              STBY(L)=MIN(CDMAXV,STBY(L))  
            ENDIF  
          ENDIF  
        ENDDO  
!$OMP END PARALLEL

      ELSEIF(N.EQ.-2)THEN
        DO L=2,LA  
          LS=LSC(L)  
          ZBRATU=0.5*(DXP(LWEST(L))*ZBR(LWEST(L))+DXP(L)*ZBR(L))*DXIU(L)  
          ZBRATV=0.5*(DYP(LS )*ZBR(LS )+DYP(L)*ZBR(L))*DYIV(L)  
          UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L)+1.E-12 )  
          VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1)+1.E-12 )  
          CDMAXU=CDLIMIT*STBXO(L)*H1U(L)/( DELT*UMAGTMP )  
          CDMAXV=CDLIMIT*STBYO(L)*H1V(L)/( DELT*VMAGTMP )  
          HURTMP=MAX(ZBRATU,H1U(L))  
          HVRTMP=MAX(ZBRATV,H1V(L))  
          DZHUDZBR=1.+0.5*DZC(1)*HURTMP/ZBRATU  
          DZHVDZBR=1.+0.5*DZC(1)*HVRTMP/ZBRATV  
!  
          STBX(L)=AVCON*STBXO(L)*.16/((LOG(DZHUDZBR))**2)  
          STBY(L)=AVCON*STBYO(L)*.16/((LOG(DZHVDZBR))**2)  
          STBX(L)=MIN(CDMAXU,STBX(L))  
          STBY(L)=MIN(CDMAXV,STBY(L)) 
        ENDDO  
      ENDIF  
!  
! **  END ROUGH DRAG FORMULATION  
!  

!
      IF(ISVEG.GE.1)THEN  
        DO K=1,KC  
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              M=MVEGL(L)
	      IF(M>90)M=M-90 !SCJ vegetation below MHK device, which is OK, it  acccounts for the structure
              FXVEG(L,K)=0.  
              FYVEG(L,K)=0.  
!  
! *** DSLLC BEGIN BLOCK  
!  
              IF(M.NE.MVEGOW.AND.M.NE.0)THEN
                ! *** M=0 FOR OPEN WATER  
                LW=LWEST(L)  
                LE=LEAST(L)  
                LS=LSC(L)  
                LN=LNC(L)  
                LNW=LNWC(L)  
                LSE=LSEC(L)  
                MW=MVEGL(LW)
	        IF(MW>90)MW=MW-90 !SCJ, vegetation below MHK device, which is OK, it  acccounts for the structure
                MS=MVEGL(LS)  
	        IF(MS>90)MS=MS-90 !SCJ, vegetation below MHK device, which is OK, it  acccounts for the structure
                VTMPATU=0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))  
                UTMPATV=0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))  
                UMAGTMP=SQRT( U(L,K)*U(L,K)+VTMPATU*VTMPATU +1.E-12 )  
                VMAGTMP=SQRT( UTMPATV*UTMPATV+V(L,K)*V(L,K) +1.E-12 )  
                UMAGTMP=MAX(UMAGTMP,UVEGSCL)  
                VMAGTMP=MAX(VMAGTMP,UVEGSCL)  
                CDMAXU=CDLIMIT*STBXO(L)*H1U(L)/( DELT*UMAGTMP )  
                CDMAXV=CDLIMIT*STBYO(L)*H1V(L)/( DELT*VMAGTMP )  
                IF(N.EQ.-2)THEN  
                  VTMPATU=0.25*(V1(L,K)+V1(LW,K)+V1(LN,K)+V1(LNW,K))  
                  UTMPATV=0.25*(U1(L,K)+U1(LE,K)+U1(LS,K)+U1(LSE,K))  
                  UMAGTMP=SQRT( U1(L,K)*U1(L,K)+VTMPATU*VTMPATU+1.E-12 )  
                  VMAGTMP=SQRT( UTMPATV*UTMPATV+V1(L,K)*V1(L,K)+1.E-12 )  
                ENDIF  
!JH             CPVEGU=0.5 ! CHANGED DEFINITION
                CPVEGU=1.0
                IF(ISVEGL.EQ.1) CPVEGU=CPVEGU + 10.E-6/(  
     &              (BPVEG(MW)+BPVEG(M))*UMAGTMP )  
                IF(CPVEGU.GT.1.0)THEN  
!  
!            CALCULATE R FOR LAMINAR FLOW  
!  
                  CPVEGU=CPVEGU-0.5  
                  RVEGUM=0.  
                ENDIF  
                CPVEGU=SCVEG(M)*CPVEGU  
!JH             CPVEGV=0.5 ! CHANGED DEFINITION
                CPVEGV=1.0
                IF(ISVEGL.EQ.1) CPVEGV=CPVEGV + 10.E-6/(  
     &              (BPVEG(MS)+BPVEG(M))*VMAGTMP )  
                IF(CPVEGV.GT.1.0)THEN  
!  
!            CALCULATE R FOR LAMINAR FLOW  
!  
                  CPVEGV=CPVEGV-0.5  
                  RVEGVM=0.  
                ENDIF  
                CPVEGV=SCVEG(M)*CPVEGV  
                FRACLAY=FLOAT(K)/FLOAT(KC)  
                FHLAYC=FRACLAY*HP(L)  
                FHLAYW=FRACLAY*HP(LWEST(L))  
                FHLAYS=FRACLAY*HP(LS)  
                HVGTC=HP(L)  
                HVGTW=HP(LWEST(L))  
                HVGTS=HP(LS)  
                IF(HPVEG(M).LT.FHLAYC) HVGTC=0.0  
                IF(HPVEG(MW).LT.FHLAYW) HVGTW=0.0  
                IF(HPVEG(MS).LT.FHLAYS) HVGTS=0.0  
               
               FXVEG(L,K)=0.25*CPVEGU*(DXP(L)*(BDLPSQ(M)*HVGTC/PVEGZ(M))
     &              +DXP(LW)*(BDLPSQ(MW)*HVGTW/PVEGZ(MW)) )*DXIU(L) 
      
               FYVEG(L,K)=0.25*CPVEGV*(DYP(L)*(BDLPSQ(M)*HVGTC/PVEGZ(M))
     &              +DYP(LS)*(BDLPSQ(MS)*HVGTS/PVEGZ(MS)) )*DYIV(L) 
                FXVEG(L,K)=MIN(FXVEG(L,K),CDMAXU)  
                FYVEG(L,K)=MIN(FYVEG(L,K),CDMAXU)
               ENDIF  
!  
! *** DSLLC END BLOCK  
!  
            ENDIF  
          ENDDO  
        ENDDO  
      ENDIF  
  300 CONTINUE  
!  
! ** SUBGRID SCALE CHANNEL FRICTION  
!  
      IF(MDCHH.GE.1)THEN  
        DO NMD=1,MDCHH  
          LHOST=LMDCHH(NMD)  
          LCHNU=LMDCHU(NMD)  
          LCHNV=LMDCHV(NMD)  
          MH=MVEGL(LHOST)  
!  
!         X-DIRECTION CHANNEL  
!  
          IF(MDCHTYP(NMD).EQ.1)THEN  
            MU=0  
            IF(ISVEG.GE.1) MU=MVEGL(LCHNU)  
            WCHAN=DXP(LCHNU)  
            RLCHN=0.5*DYP(LCHNU)+CHANLEN(NMD)  
            HCHAN=0.5*DYP(LCHNU)*H1P(LCHNU)+CHANLEN(NMD)*H1P(LHOST)  
            HCHAN=HCHAN/RLCHN  
            ZBRATU=0.5*DYP(LCHNU)*ZBR(LCHNU)+CHANLEN(NMD)*ZBR(LHOST)  
            ZBRATU=ZBRATU/RLCHN  
            HURTMP=MAX(ZBRATU,HCHAN)  
            HUDZBR=HURTMP/ZBRATU  
            IF(HUDZBR.LT.7.5) HUDZBR=7.5  
            STBXCH=0.16/( (LOG( HUDZBR ) -1.)**2)  
            CDMAXU=HCHAN*HCHAN*WCHAN/( DELT*(QCHANU(NMD)+1.E-12) )  
            STBXCH=MAX(STBXCH,CDMAXU)  
            STBXCH=MAX(STBXCH,0.1)  
            FXVEGCH=0.0  
            IF(MU.GT.0) FXVEGCH=  
     &          0.5*(0.5*DYP(LCHNU)*(BDLPSQ(MU)*H1P(LCHNU)/PVEGZ(MU))  
     &          +CHANLEN(NMD)*(BDLPSQ(MH)*H1P(LHOST)/PVEGZ(MH)) )/RLCHN  
            CHANFRIC(NMD)=FXVEGCH+STBXCH  
          ENDIF  
!  
!         Y-DIRECTION CHANNEL  
!  
          IF(MDCHTYP(NMD).EQ.2)THEN  
            MV=0  
            IF(ISVEG.GE.1) MV=MVEGL(LCHNV)  
            WCHAN=DYP(LCHNV)  
            RLCHN=0.5*DXP(LCHNV)+CHANLEN(NMD)  
            HCHAN=0.5*DXP(LCHNV)*H1P(LCHNV)+CHANLEN(NMD)*H1P(LHOST)  
            HCHAN=HCHAN/RLCHN  
            ZBRATV=0.5*DXP(LCHNV)*ZBR(LCHNV)+CHANLEN(NMD)*ZBR(LHOST)  
            ZBRATV=ZBRATV/RLCHN  
            HVRTMP=MAX(ZBRATV,HCHAN)  
            HVDZBR=HVRTMP/ZBRATV  
            IF(HVDZBR.LT.7.5) HVDZBR=7.5  
            STBYCH=0.16/( (LOG( HVDZBR ) -1.)**2)  
            CDMAXV=HCHAN*HCHAN*WCHAN/( DELT*(QCHANV(NMD)+1.E-12) )  
            STBYCH=MAX(STBYCH,CDMAXV)  
            STBYCH=MAX(STBYCH,0.1)  
            FYVEGCH=0.0  
            IF(MV.GT.0) FYVEGCH=  
     &          0.5*(0.5*DXP(LCHNV)*(BDLPSQ(MV)*H1P(LCHNV)/PVEGZ(MV))  
     &          +CHANLEN(NMD)*(BDLPSQ(MH)*H1P(LHOST)/PVEGZ(MH)) )/RLCHN  
            CHANFRIC(NMD)=FYVEGCH+STBYCH  
          ENDIF  
        ENDDO  
      ENDIF  
      IF(ISVEG.GE.2.AND.KC.GT.1)THEN  
        DO L=2,LA  
          M=MVEGL(L)  
          MW=MVEGL(LWEST(L))  
          MS=MVEGL(LSC(L))  
          WRITE(1,1122)N,IL(L),JL(L),MVEGL(L),PVEGZ(M),PVEGZ(MS),  
     &        PVEGZ(MW),STBX(L),STBY(L)  
          WRITE(1,1123)(FXVEG(L,K),K=1,KC)  
          WRITE(1,1123)(FYVEG(L,K),K=1,KC)  
        ENDDO  
      ENDIF  
      IF(ISVEG.GE.2) CLOSE(1)  
 1122 FORMAT(4I5,5E12.4)  
 1123 FORMAT(15X,10E12.4)  
      GOTO 1948  
!  
! **  ENTER HERE FOR WAVE-CURRENT BOUNDARY LAYER  
!  
 1947 CONTINUE  
      IF(JSTBXY.EQ.0)THEN  
        DO L=2,LA  
          STBXO(L)=STBX(L)  
          STBYO(L)=STBY(L)  
        ENDDO  
        N=0  
        JSTBXY=1  
        IF(ISDZBR.GE.1)THEN  
          OPEN(1,FILE='ZBREMX.OUT',STATUS='UNKNOWN')  
          CLOSE(1,STATUS='DELETE')  
        ENDIF  
      ENDIF  
      IF(ISDZBR.EQ.N.AND.DEBUG)THEN  
        OPEN(1,FILE='CDDIAG.OUT',STATUS='UNKNOWN')  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE='CDDIAG.OUT',STATUS='UNKNOWN')  
      ENDIF  
      NTMP=MAX(N,1)  
      IF(NTMP.LT.NTSWV)THEN  
        TMPVALW=FLOAT(NTMP)/FLOAT(NTSWV)  
        WVFACT=0.5-0.5*COS(PI*TMPVALW)  
      ELSE  
        WVFACT=1.0  
      ENDIF  
!  
! *** DSLLC BEGIN BLOCK  
!  
      DO L=2,LA  
        IF(UWVSQ(L).GT.1.E-6 .AND. LMASKDRY(L))THEN  
          QQWCTMP=SQRT( QQWV2(L)*QQWV2(L)+QQ(L,0)*QQ(L,0) )  
          TWCTMP=QQWCTMP/CTURB2  
!  
!        CORZBR=1.+1.2*TAUTMP/(1.+0.2*TAUTMP)  
!  
          AEXTMP=WVWHA(L)/SINH(WVKHP(L))  
          ZBRE(L)=ZBR(L)  
          IF(QQ(L,0).GT.0.)THEN  
            TMPVAL=UWVSQ(L)*SQRT( AEXTMP/(30.*ZBR(L)) )  
            USTARC=SQRT(QQ(L,0)/CTURB2)  
            TMPVAL=TMPVAL/USTARC  
            ZBRE(L)=ZBR(L)*(1.+0.19*TMPVAL)  
          ENDIF  
          CDRGTMP=(30.*ZBRE(L)/AEXTMP)**0.2  
          CDRGTMP=5.57*CDRGTMP-6.13  
          CDRGTMP=EXP(CDRGTMP)  
          CDRGTMP=MIN(CDRGTMP,0.22)  
          TAUTMP=0.5*CDRGTMP*UWVSQ(L)  
          QQWV2(L)=CTURB2*TAUTMP*WVFACT  
          QQWC(L)=SQRT( QQWV2(L)*QQWV2(L)+QQ(L,0)*QQ(L,0) )  
          IF(ISTRAN(7).GT.0)THEN  
            TWCTMP=QQWC(L)/CTURB2  
            TAUBTMP=QQWV1(L)/CTURB2  
            TAUE=TWCTMP/TAUN(NSED+1)  
            RIPAMP=0.  
            RIPSTP=0.  
            IF(TAUBTMP.GT.TAUN(NSED+1).AND.TAUBTMP.LE.TAUD(NSED+1))THEN  
              RIPAMP=0.22/(TAUE**0.16)  
              RIPSTP=0.16/(TAUE**0.04)  
            ENDIF  
            IF(TAUBTMP.GT.TAUD(NSED+1))THEN  
              RIPAMP=0.78/(TAUE**1.5)  
              RIPSTP=0.41/TAUE  
            ENDIF  
            RIPAMP=RIPAMP*WVWHA(L)/SINH(WVKHP(L))  
            TMPVAL=0.  
            IF(RIPAMP.GT.0.) TMPVAL=LOG(RIPAMP/ZBRE(L))-1.  
            TMPVAL=MAX(TMPVAL,0.)  
            RIPFAC=1.+3.125*TMPVAL*TMPVAL*RIPSTP  
            QQWV3(L)=RIPFAC*QQWV2(L)  
            QQWCR(L)=SQRT( QQWV3(L)*QQWV3(L)+QQ(L,0)*QQ(L,0) )  
          ELSE  
            QQWCR(L)=QQ(L,0)  
          ENDIF  
        ELSE  
          QQWV2(L)=QQLMIN  
          QQWC(L)=QQ(L,0)  
          QQWCR(L)=QQ(L,0)  
        ENDIF  
      ENDDO  
!  
! *** DSLLC END BLOCK  
!  
      ZBRMAX=-(1.E+12)*ZBRADJ  
      ZBRMIN=(1.E+12)*ZBRADJ  
      CDRGMAX=-1.E+12  
      CDRGMIN=1.E+12  
      IF(ISWAVE.EQ.1.OR.ISWAVE.EQ.2)WVDTMP=0.4/(WVFRQ*CTURB3)  
      RKZTURB=0.4/CTURB3  
      DO L=2,LA  
        IF(LMASKDRY(L))THEN
          LS=LSC(L)  
          LN=LNC(L)  
          UTMP=0.5*STCUV(L)*(U(LEAST(L),1)+U(L,1))+1.E-12  
          VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))  
          CURANG=ATAN2(VTMP,UTMP)  
          COSWC=COS(CURANG-WACCWE(L))  
          UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L)+1.E-12 )  
          VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1)+1.E-12 )  
          CDMAXU=STBXO(L)*H1U(L)/( 4.*DELT*UMAGTMP )  
          CDMAXV=STBYO(L)*H1V(L)/( 4.*DELT*VMAGTMP )  
          CDTMPU=-1.  
          CDTMPV=-1.  
          QWCTMPU=0.5*( QQWV2(L)+QQWV2(LEAST(L)) )  
          QWCTMPV=0.5*( QQWV2(L)+QQWV2(LS ) )  
          IF(ISWCBL.EQ.2)THEN  
            QWCTMPU=0.5*( QQWC(L)+QQWC(LEAST(L)) )  
            QWCTMPV=0.5*( QQWC(L)+QQWC(LS ) )  
          ENDIF 
          IF(ISWAVE.EQ.3)THEN  
            IF(WVFRQL(L).GT.1E-6)THEN
              WVDTMP=0.4/(WVFRQL(L)*CTURB3)
            ELSE
              WVDTMP=0.
            ENDIF
          ENDIF  
          WVDELU=WVDTMP*SQRT(QWCTMPU)  
          WVDELV=WVDTMP*SQRT(QWCTMPV)  
          QWCTMPU=0.5*( QQWCR(L)+QQWCR(LEAST(L)) )  
          QWCTMPV=0.5*( QQWCR(L)+QQWCR(LS ) )  
          QWCTMPU=SQRT(QWCTMPU)  
          QWCTMPV=SQRT(QWCTMPV)  
          QCTMPU=0.5*( QQ(L,0)+QQ(LEAST(L),0) )  
          QCTMPV=0.5*( QQ(L,0)+QQ(LS ,0) )  
          QWDQCU=QWCTMPU/SQRT(QCTMPU)  
          QWDQCV=QWCTMPV/SQRT(QCTMPV)  
          HZREFU=DZC(1)*H1U(L)  
          HZREFV=DZC(1)*H1V(L)  
          ZBREU=0.5*(ZBRE(L)+ZBRE(LEAST(L)))  
          ZBREV=0.5*(ZBRE(L)+ZBRE(LS ))  
          ZDHZRU=ZBREU/HZREFU  
          ZDHZRV=ZBREV/HZREFV  
          HZRUDZ=1./ZDHZRU  
          HZRVDZ=1./ZDHZRV  
          DWUD2Z=0.5*WVDELU/ZBREU  
          DWVD2Z=0.5*WVDELV/ZBREV  
          DWUDZ=2.*DWUD2Z  
          DWVDZ=2.*DWVD2Z  
          DWUDHR=WVDELU/HZREFU  
          DWVDHR=WVDELV/HZREFV  
          CDTMPUX=RKZTURB*QWCTMPU  
          CDTMPVY=RKZTURB*QWCTMPV  
          JWCBLU=0  
          JWCBLV=0  
          IF( HZRUDZ.LE.DWUD2Z)THEN  
            CDTMPU=CDTMPUX/( (1.+ZDHZRU)*LOG(1.+HZRUDZ)-1. )  
            JWCBLU=1  
          ENDIF  
          IF( HZRVDZ.LE.DWVD2Z)THEN  
            CDTMPV=CDTMPVY/( (1.+ZDHZRV)*LOG(1.+HZRVDZ)-1. )  
            JWCBLV=1  
          ENDIF  
          IF( HZRUDZ.GT.DWUD2Z.AND.HZRUDZ.LE.DWUDZ)THEN  
            BOTTMP=(1.+ZDHZRU)*LOG(1.+DWUD2Z)-0.5*DWUDHR  
     &        +0.5*HZRUDZ*(1.-0.5*DWUDHR)*(1.-0.5*DWUDHR)/(1.+DWUD2Z)  
          CDTMPU=CDTMPUX/BOTTMP  
          JWCBLU=2  
        ENDIF  
        IF( HZRVDZ.GT.DWVD2Z.AND.HZRVDZ.LE.DWVDZ)THEN  
          BOTTMP=(1.+ZDHZRV)*LOG(1.+DWVD2Z)-0.5*DWVDHR  
     &        +0.5*HZRVDZ*(1.-0.5*DWVDHR)*(1.-0.5*DWVDHR)/(1.+DWVD2Z)  
          CDTMPV=CDTMPVY/BOTTMP  
          JWCBLV=2  
        ENDIF  
        IF( HZRUDZ.GT.DWUDZ)THEN  
          BOTTMP=QWDQCU*( (1.+ZDHZRU)*(LOG(1.+HZRUDZ)-LOG(1.+DWUDZ))  
     &        +DWUDHR-1. )  
          BOTTMP=BOTTMP+(1.+ZDHZRU)*LOG(1.+DWUD2Z)  
     &        +DWUD2Z*(1.-1.25*DWUDHR-ZDHZRU)/(1.+DWUD2Z)  
          CDTMPU=CDTMPUX/BOTTMP  
          JWCBLU=3  
        ENDIF  
        IF( HZRVDZ.GT.DWVDZ)THEN  
          BOTTMP=QWDQCV*( (1.+ZDHZRV)*(LOG(1.+HZRVDZ)-LOG(1.+DWVDZ))  
     &        +DWVDHR-1. )  
          BOTTMP=BOTTMP+(1.+ZDHZRV)*LOG(1.+DWVD2Z)  
     &        +DWVD2Z*(1.-1.25*DWVDHR-ZDHZRV)/(1.+DWVD2Z)  
            CDTMPV=CDTMPVY/BOTTMP  
            JWCBLV=3  
          ENDIF  
          CDTMPU=CDTMPU/UMAGTMP  
          CDTMPV=CDTMPV/VMAGTMP
          IF(DEBUG)THEN  
            IF(ISDZBR.EQ.N)THEN  
              WRITE(1,1779) IL(L),JL(L),JWCBLU,JWCBLV  
              WRITE(1,1780)  
              WRITE(1,1781) ZBREU,WVDELU,HZREFU,CDTMPU,CDMAXU  
              WRITE(1,1782)  
              WRITE(1,1781) ZBREV,WVDELV,HZREFV,CDTMPV,CDMAXV  
            ENDIF  
          ENDIF
          IF(CDTMPU.LE.0.) CDTMPU=CDMAXU  
          IF(CDTMPV.LE.0.) CDTMPV=CDMAXV  
          STBX(L)=AVCON*STBXO(L)*CDTMPU  
          STBY(L)=AVCON*STBYO(L)*CDTMPV  
          STBX(L)=MIN(CDMAXU,STBX(L),0.11)  
          STBY(L)=MIN(CDMAXV,STBY(L),0.11) 
        ENDIF
      ENDDO
      IF(DEBUG)THEN  
        IF(ISDZBR.EQ.N) CLOSE(1)  
        IF(ISDZBR.GE.1)THEN  
          DO L=2,LA  
            IF(ZBRE(L).GT.ZBRMAX)THEN  
              ZBRMAX=ZBRE(L)  
              LZBMAX=L  
            ENDIF  
            IF(ZBRE(L).LT.ZBRMIN)THEN  
              ZBRMIN=ZBRE(L)  
              LZBMIN=L  
            ENDIF  
            IF(STBX(L).GT.CDRGMAX)THEN  
              CDRGMAX=STBX(L)  
              LCDMAX=L  
            ENDIF  
            IF(STBX(L).LT.CDRGMIN)THEN  
              CDRGMIN=STBX(L)  
              LCDMIN=L  
            ENDIF  
            IF(STBY(L).GT.CDRGMAX)THEN  
              CDRGMAX=STBY(L)  
              LCDMAX=L  
            ENDIF  
            IF(STBY(L).LT.CDRGMIN)THEN  
              CDRGMIN=STBY(L)  
              LCDMIN=L  
            ENDIF  
          ENDDO  
          OPEN(1,FILE='ZBREMX.OUT',STATUS='UNKNOWN',POSITION='APPEND')  
          HOTLYMX=DZC(1)*H1P(LZBMAX)  
          HOTLYMN=DZC(1)*H1P(LZBMIN)  
          WRITE(1,1739)N,IL(LZBMAX),JL(LZBMAX),ZBRMAX,HOTLYMX  
          WRITE(1,1749)N,IL(LZBMIN),JL(LZBMIN),ZBRMIN,HOTLYMN  
          WRITE(1,1759)N,IL(LCDMAX),JL(LCDMAX),CDRGMAX,STBX(LCDMAX),  
     &      STBY(LCDMAX)  
          WRITE(1,1769)N,IL(LCDMIN),JL(LCDMIN),CDRGMIN,STBX(LCDMIN),  
     &      STBY(LCDMIN)  
          CLOSE(1)  
        ENDIF  
      ENDIF
 1948 CONTINUE  
 1717 FORMAT(' N,I,J = ',I10,2I5,'   CDTOTU,CDMAXU = ',2F15.10)  
 1718 FORMAT(' N,I,J = ',I10,2I5,'   CDTOTV,CDMAXV = ',2F15.10)  
 1727 FORMAT(' N,I,J = ',I10,2I5,'   LAM CDTOTU,CDMAXU = ',2F15.10)  
 1728 FORMAT(' N,I,J = ',I10,2I5,'   LAM CDTOTV,CDMAXV = ',2F15.10)  
 1719 FORMAT(' N = ',I10,'  CDTOTUM,CDTOTVM = ',2F15.10)  
 1729 FORMAT(' N = ',I10,'  CDMAXUM,CDMAXVM = ',2F15.10)  
 1739 FORMAT(' N,I,J = ',I10,2I5,'  ZBRMAX,HBTLYMX = ',2E14.6)  
 1749 FORMAT(' N,I,J = ',I10,2I5,'  ZBRMIN,HBTLYMN = ',2E14.6)  
 1759 FORMAT(' N,I,J = ',I10,2I5,'  CDRGMAX,STBX,STBY = ',3E14.6)  
 1769 FORMAT(' N,I,J = ',I10,2I5,'  CDRGMIN,STBX,STBY = ',3E14.6)  
 1779 FORMAT(' I, J, JWCBLU, JWCBLV = ',4I8)  
 1780 FORMAT('    ZBREU        WVDELU        HZREFU        CDTMPU    ',  
     &    1X,'  CDMAXU')  
 1781 FORMAT(5E12.4)  
 1782 FORMAT('    ZBREV        WVDELV        HZREFV        CDTMPV    ',  
     &    1X,'  CDMAXV')  
      RETURN  
      END  

