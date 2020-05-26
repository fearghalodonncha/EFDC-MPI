      SUBROUTINE AINIT  
C  
C CHANGE RECORD  
C  ADDED TRANSPORT BYPASS MASK, IMASKDRY FOR DRY CELLS  
C  ADDED TRANSPORT BYPASS MASK, LMASKDRY FOR DRY CELLS  
C  MODIFIED DEFINITION OF CHANLEN IN INITIALIZATION RATHER THAN  
C  IN SUBS CALTBXY AND CALPUV2C AND CALPUV9C  
C
C  ALL ZEROING OF ARRAYS MOVED TO ZERO
C  
      USE GLOBAL  
	IMPLICIT NONE
	INTEGER::L,I,J,LS,NT,LCHNV,IVAL,NS,K,NMD,LHOST,LCHNU,NV,NX
	INTEGER::NTMPC,NTMPN
C  
C **  INITIALIZE ARRAYS  
C  
      ZBR(1)=ZBRADJ  
      ZBRE(1)=ZBRADJ  
      HMP(1)=HMIN  
      HMU(1)=HMIN  
      HMV(1)=HMIN  
      HWQ(1)=HMIN  
      H2WQ(1)=HMIN  
      DXP(1)=DX  
      DYP(1)=DY  
      DXU(1)=DX  
      DYU(1)=DY  
      DXV(1)=DX  
      DYV(1)=DY  
      DXYP(1)=DX*DY  
      MVEGL(1)=0
      BELV(1)=BELV(2)  
      ZBR(LC)=ZBRADJ  
      ZBRE(LC)=ZBRADJ  
      HMP(LC)=HMIN  
      HMU(LC)=HMIN  
      HMV(LC)=HMIN  
      HWQ(LC)=HMIN  
      H2WQ(LC)=HMIN  
      DXP(LC)=DX  
      DYP(LC)=DY  
      DXU(LC)=DX  
      DYU(LC)=DY  
      DXV(LC)=DX  
      DYV(LC)=DY  
      DXYP(LC)=DX*DY  
      MVEGL(LC)=0  
      BELV(LC)=BELV(LA)  
      IF(ISGWIE.EQ.0) DAGWZ=0.0
      DO L=2,LA  
        I=IL(L)  
        J=JL(L)
        KBT(L)=1
        BELAGW(L)=BELV(L)-DAGWZ  
        ZBRE(L)=ZBR(L)  
        DLON(L)=CDLON1+(CDLON2*FLOAT(I)+CDLON3)/60.0
        DLAT(L)=CDLAT1+(CDLAT2*FLOAT(J)+CDLAT3)/60.0
        CUE(L)=1.0
        CVN(L)=1.0 
      ENDDO  
      HMU(1)=HMU(2)    ! *** PMC
      HMV(1)=HMV(2)    ! *** PMC
      HMU(LC)=HMU(LA)  ! *** PMC
      HMV(LC)=HMV(LA)  ! *** PMC
      DO L=1,LC  
        CC(L)=1.0
        CCC(L)=1.0
        P(L)=G*(HMP(L)+BELV(L))  
        P1(L)=G*(HMP(L)+BELV(L))  
        HP(L)=HMP(L)+PDGINIT  
        HU(L)=HMU(L)+PDGINIT  
        HV(L)=HMV(L)+PDGINIT  
        HPI(L)=1.0/HP(L)  
        HUI(L)=1.0/HU(L)  
        HVI(L)=1.0/HV(L)  
        HWQ(L)=HMP(L)+PDGINIT  
        H1P(L)=HMP(L)+PDGINIT  
        H2P(L)=HMP(L)+PDGINIT  
        H1U(L)=HMU(L)+PDGINIT  
        H1V(L)=HMV(L)+PDGINIT  
        H1UI(L)=1./H1U(L)  
        H1VI(L)=1./H1V(L)  
        H2WQ(L)=HMP(L) !+PDGINIT  
        SCB(L)=1.0
        SPB(L)=1.0
        SUB(L)=1.0
        SVB(L)=1.0
        SWB(L)=1.0
        STCUV(L)=1.0
        STCAP(L)=1.0
        STBX(L)=1.0
        STBY(L)=1.0
        SAAX(L)=1.0
        SAAY(L)=1.0
        SNLPX(L)=1.0
        SNLPY(L)=1.0
        SCAX(L)=1.0
        SCAY(L)=1.0
        SBX(L)=1.0
        SBY(L)=1.0
        SDX(L)=1.0
        SDY(L)=1.0
        LMASKDRY(L)=.TRUE.  
      ENDDO  
C
C *** DSLLC BEGIN BLOCK
      ! *** OPEN WATER DEFAULT SETTINGS
      NV=0
      PVEGX(NV)=1.0
      PVEGY(NV)=1.0
      PVEGZ(NV)=1.0
C *** DSLLC END BLOCK
C
      DO NT=1,NTOX  
        DO K=1,KB  
          DO L=1,LC  
            TOXB(L,K,NT)=TOXBINIT(L,K,NT)  
            TOXB1(L,K,NT)=TOXBINIT(L,K,NT)  
          ENDDO  
        ENDDO  
      ENDDO
      IF(IWRSP(1).LT.98)THEN
        DO NS=1,NSED  
          DO K=1,KB  
            DO L=1,LC  
              SEDB(L,K,NS)=SEDBINIT(L,K,NS)  
              SEDB1(L,K,NS)=SEDBINIT(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDDO
      ENDIF
      DO NS=1,NSND  
        NX=NS+NSED  
        DO K=1,KB  
          DO L=1,LC  
            SNDB(L,K,NS)=SNDBINIT(L,K,NS)  
            SNDB1(L,K,NS)=SNDBINIT(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  

C      IF(IS1DCHAN.EQ.1)THEN  
C        DO L=1,LC  
C          FADYP(L)=1.  
C          FADYP1(L)=1.  
C          FADYP2(L)=1.  
C          WPDYP(L)=1.  
C          WPDYP1(L)=1.  
C          FADXP(L)=1.  
C          FADXP1(L)=1.  
C          FADXP2(L)=1.  
C          WPDXP(L)=1.  
C          WPDXP1(L)=1.  
C          FADYU(L)=1.  
C          FADYU1(L)=1.  
C          WPDYU(L)=1.  
C          WPDYU1(L)=1.  
C          FADXV(L)=1.  
C          FADXV1(L)=1.  
C          WPDXV(L)=1.  
C          WPDXV1(L)=1.  
C          DADH(L)=1.  
C          DADH1(L)=1.  
C          SRFXP(L)=0.  
C          SRFYP(L)=0.  C
C          SRFXP1(L)=0.  
C          SRFYP1(L)=0.  
C          SRFXV(L)=0.  
C          SRFYU(L)=0.  
C          SRFXV1(L)=0.  
C          SRFYU1(L)=0.  
C        ENDDO  
C      ENDIF  
      DO L=1,NLRPD  
        NLRPDL(L)=1  
      ENDDO  
      DO K=1,KS  
        DO L=1,LC  
          AV(L,K)=AVO  
          AVVI(L,K)=1./AVO  
          AVUI(L,K)=1./AVO  
          AB(L,K)=ABO  
          QQL(L,K)=QQLMIN  
          QQL1(L,K)=QQLMIN  
          QQL2(L,K)=QQLMIN  
          DML(L,K)=DMLMIN  
C *** ALL ZEROING OF ARRAYS MOVED TO ZERO
        ENDDO  
      ENDDO  
!      DO K=1,KC  
!        DO L=1,LC  
          AH(1:LC,1:KC)=AHO  
          AHU(1:LC,1:KC)=AHO  
          AHULPF(1:LC,1:KC)=AHO  
          AHV(1:LC,1:KC)=AHO  
          AHVLPF(1:LC,1:KC)=AHO  
          AHC(1:LC,1:KC)=AHO  
          AQ(1:LC,1:KC)=AVO  
C *** ALL ZEROING OF ARRAYS MOVED TO ZERO
          CTURBB1(1:LC,1:KC)=CTURB  
          CTURBB2(1:LC,1:KC)=CTURB2B  
C *** TEMPERATURE INITIATION
          TEM(1:LC,1:KC)=TEMO    !This shouldn't be here, it is already initialized in INPUT.f...SCJ...where?
          TEM1(1:LC,1:KC)=TEM(1:LC,1:KC)   
!        ENDDO  
!      ENDDO  
      NTMPC=MAX(NSED,1)  
      DO NS=1,NTMPC  
        DO K=1,KC  
          DO L=1,LC  
            SED(L,K,NS)=SEDO(NS)  
            SED1(L,K,NS)=SEDO(NS)  
C *** ALL ZEROING OF ARRAYS MOVED TO ZERO
          ENDDO  
        ENDDO  
      ENDDO  
      NTMPN=MAX(NSND,1)  
      DO NX=1,NTMPN  
        NS=NX+NTMPC  
        DO K=1,KC  
          DO L=1,LC  
            SND(L,K,NX)=SEDO(NS)  
            SND1(L,K,NX)=SEDO(NS)  
C *** ALL ZEROING OF ARRAYS MOVED TO ZERO
          ENDDO  
        ENDDO  
      ENDDO  
      DO NT=1,NTOX  
        DO K=1,KC  
          DO L=1,LC  
            TOX(L,K,NT)=TOXINTW(NT)  
            TOX1(L,K,NT)=TOXINTW(NT)  
C *** ALL ZEROING OF ARRAYS MOVED TO ZERO
          ENDDO  
        ENDDO  
      ENDDO  
      DO K=0,KC  
        DO L=1,LC  
C *** ALL ZEROING OF ARRAYS MOVED TO ZERO
          QQ(L,K)=QQMIN  
          QQ1(L,K)=QQMIN  
          QQ2(L,K)=QQMIN
          QQSQR(L,K)=SQRT(QQMIN)
        ENDDO  
      ENDDO  
      IF(MDCHH.GE.1)THEN  
        DO NMD=1,MDCHH  
          LHOST=LMDCHH(NMD)  
          LCHNU=LMDCHU(NMD)  
          LCHNV=LMDCHV(NMD)  
C  
C         SET HOST DRYING DEPTH  
C  
          IF(PMDCH(NMD).LT.0.0) PMDCH(NMD)=HWET  
C  
C         X-DIRECTION CHANNEL  
C  
          IF(MDCHTYP(NMD).EQ.1)THEN  
            IF(CHANLEN(NMD).LT.0.0)THEN  
              CHANLEN(NMD)=0.25*DYP(LHOST)  
            ELSE  
              CHANLEN(NMD)=CHANLEN(NMD)-0.5*DYP(LCHNU)  
            ENDIF  
          ENDIF  
C  
C         Y-DIRECTION CHANNEL  
C  
          IF(MDCHTYP(NMD).EQ.2)THEN  
            IF(CHANLEN(NMD).LT.0.0)THEN  
              CHANLEN(NMD)=0.25*DXP(LHOST)  
            ELSE  
              CHANLEN(NMD)=CHANLEN(NMD)-0.5*DXP(LCHNV)  
            ENDIF  
          ENDIF  
        ENDDO  
      ENDIF  
C  
C ** INITIALIZE ORGANIC CARBON VARIABLES OF SEDIMENT-TOXICS  
C  
      IVAL=0  
      DO NT=1,NTOX  
        IF(ISTOC(NT).GT.0)IVAL=1  
      ENDDO  
      IF(IVAL.EQ.0)THEN  
        DO NS=1,NSED+NSND  
          DO K=1,KB  
            DO L=1,LC  
              STFPOCB(L,K,NS)=1.0  
            ENDDO  
          ENDDO  
        ENDDO  
        DO NS=1,NSED+NSND  
          DO K=1,KC  
            DO L=1,LC  
              STFPOCW(L,K,NS)=1.0  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  
      IF(IVAL.EQ.1)THEN  
        IF(ISTDOCB.EQ.0)THEN  
          DO K=1,KB  
            DO L=1,LC  
              STDOCB(L,K)=STDOCBC  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISTPOCB.EQ.0)THEN  
          DO K=1,KB  
            DO L=1,LC  
              STPOCB(L,K)=STPOCBC  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISTPOCB.EQ.2)THEN  
          DO NS=1,NSED+NSND  
            DO K=1,KB  
              DO L=1,LC  
                STFPOCB(L,K,NS)=FPOCBST(NS,1)  
              ENDDO  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISTDOCW.EQ.0)THEN  
          DO K=1,KC  
            DO L=1,LC  
              STDOCW(L,K)=STDOCWC  
            ENDDO  
          ENDDO  
        ENDIF  
        IF(ISTPOCW.EQ.0)THEN  
          DO K=1,KC  
            DO L=1,LC  
              STPOCW(L,K)=STPOCWC  
            ENDDO
          ENDDO
        ENDIF
        IF(ISTPOCW.EQ.2)THEN
          DO NT=1,NTOX
            DO NS=1,NSED+NSND  
              DO K=1,KC  
                DO L=1,LC  
                  STFPOCW(L,K,NS)=FPOCWST(NS,NT)  
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C**********************************************************************C
CGR 4/16/10 INITIALIZATION OF PUVDASM LAYER CONSTRAINTS
C     CONSTRAIN ASSIMILATION TO DEPTHS ABOVE KCUVDA W/SMOOTHING
      TUVKC(1:KC)=1.0
      DO K=1,KCUVDA
        TUVKC(K)=TUVKC(K)*FLOAT(K)/FLOAT(KCUVDA)
      ENDDO  
c CSERT      
      TCSERKC(1:KC)=1.0
      DO K=1,KCSERDA
        TCSERKC(K)=TCSERKC(K)*FLOAT(K)/FLOAT(KCSERDA)
      ENDDO  
C********************* 
C***********************************************************	
CGR 6/8/10  INITIALIZE CDASM
C ***Relax CSERext. Optimized for sponge width of 5-10 cells.  No Sponge Layers set by NLCDA in C16***
C****Greg Rocheleau       Feb 10,2010  
C
CGR 5/10/2011      IF(ISACDA.EQ.3)THEN
CGR 5/10/2011        DO NSPNG=1,NLCDA    
C 6/8/10 ORIGINAL CODE WAS "SWAPPING" INDEXING OS R1, BUT WAS WORKING WELL...
C       R1(NSPNG)=TSCDA*(REAL(NSPNG)/REAL(NLCDA))**10
CGR 5/10/2011        R1(NLCDA-NSPNG+1)=TSCDA*(REAL(NSPNG)/REAL(NLCDA))**6
C	  R2(NSPNG)=1.
CGR 5/10/2011        R2(NLCDA-NSPNG+1)=TSCDA*(REAL(NSPNG)/REAL(NLCDA))**3
CGR 5/10/2011	  ENDDO
C
C        DO NSPNG=NLCDA+1,NLCDA+3    
C	  R2(NSPNG)=(0.5)**(REAL(NSPNG)-REAL(NLCDA))
C	  ENDDO
CGR 5/10/2011	ENDIF
	IF(ISACDA==3.AND.ISCDA(3)>0)R1(1:ISCDA(3))=1.0
C***********************************************************************
CGR
C**********************************************************************C	
C**********************************************************************C
      RETURN  
      END  

