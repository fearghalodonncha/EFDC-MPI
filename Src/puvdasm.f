C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE PUVDASM(ISTL_,ICALL_)
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C **  PERFORMCE WATER SURFACE ELEVATION AND HORIZONTAL VELOCITY DATA
C **  ASSIMILATION
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED     BY                 DATE APPROVED    BY
C
C----------------------------------------------------------------------C
C
C**********************************************************************C
C
      USE GLOBAL
C
      INTEGER::I,J

C
C**********************************************************************C
C

      IF(ISDYNSTP.EQ.0)THEN
        DELT=DT2
        DELTD2=DT
        IF(ISTL_.EQ.2)THEN
          DELT=DT
          DELTD2=0.5*DT
        ENDIF
        DELTI=1./DELT
      ELSE
        DELT=DTDYN
        DELTD2=0.5*DTDYN
        DELTI=1./DELT
      ENDIF
C
      NDAYA=MOD(N,NTSPTC)
      NDAYA=1+(N-NDAYA)/NTSPTC
	RVAL=-998.0
C
C**********************************************************************C
C
C **  INTERPOLATE INPUT TIME SERIES OF OBSERVED HORIZONTAL VELOCITIES
C
      DO L=1,LC
        RFBODYFXI(L)=0.0
        RFBODYFYI(L)=0.0
        FBODYFXI(L)=0.0
        FBODYFYI(L)=0.0
      ENDDO
C
      DO K=1,KC
      DO L=1,LC
	  FBODYFX(L,K)=0.0
	  FBODYFY(L,K)=0.0
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  INTERPOLATE INPUT TIME SERIES OF OBSERVED HORIZONTAL VELOCITIES
C
C----------------------------------------------------------------------C
C
      IF(ICALL_.EQ.1)THEN
      IF(NUVSER.GT.0)THEN
C
      DO K=1,KC
        USERT(K,0)=0.
        VSERT(K,0)=0.
      ENDDO
C
      DO NS=1,NUVSER
C
      IF(ISDYNSTP.EQ.0)THEN
        TIME=DT*FLOAT(N)/TCUVSER(NS)+TBEGIN*(TCON/TCUVSER(NS))
      ELSE
        TIME=TIMESEC/TCUVSER(NS)
      ENDIF
C
      M1=MUVTLAST(NS)
  100 CONTINUE
      M2=M1+1
      IF(TIME.GT.TUVSER(M2,NS))THEN
       M1=M2
       GOTO 100
      ELSE
       MUVTLAST(NS)=M1
      ENDIF      
C
      TDIFF=TUVSER(M2,NS)-TUVSER(M1,NS)
      WTM1=(TUVSER(M2,NS)-TIME)/TDIFF
      WTM2=(TIME-TUVSER(M1,NS))/TDIFF
      DO K=1,KC
        USERT(K,NS)=WTM1*USER(M1,K,NS)+WTM2*USER(M2,K,NS)
	  IF(USER(M1,K,NS).LT.RVAL)USERT(K,NS)=USER(M1,K,NS)
 	  IF(USER(M2,K,NS).LT.RVAL)USERT(K,NS)=USER(M2,K,NS)
        VSERT(K,NS)=WTM1*VSER(M1,K,NS)+WTM2*VSER(M2,K,NS)
	  IF(VSER(M1,K,NS).LT.RVAL)VSERT(K,NS)=VSER(M1,K,NS)
	  IF(VSER(M2,K,NS).LT.RVAL)VSERT(K,NS)=VSER(M2,K,NS)
	ENDDO
C
C      WRITE(7,*)TIME,NS,USER(M1,KC,NS),USER(M2,KC,NS),
C     &                  VSER(M1,KC,NS),VSER(M2,KC,NS)
C
      ENDDO
C
      ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  WATER SURFACE ELEVATION DATA ASSIMILATION
C
C     CALCULATIONS A VOLUME SOURCE OR SINK NECESSARY TO FORCE THE 
C     COMPUTE WATER SURFACE ELEVATION TOWARD THE OBSERVED WATER SURFACE
C     ELEVATION. THE VOLUME SOURCE OR SINK IS INSERTED INTO THE
C     INTERNAL AND EXTERNAL CONTINUITY EQUATIONS ON THE NEXT TIME STEP
C     THE COMPUTED SOURCE SINK FLOW IS WRITTEN TO VOLUME SOURCE
C     TIME SERIES OUTPUT AND A DAILY SUMMARY OUTPUT FILE.
C
C----------------------------------------------------------------------C
C
      IF(ICALL_.EQ.1)THEN
      IF(ISWSEDA.GT.0)THEN
C
CGR   THIS DO LOOP ADDED 6/1/2011
      DO K=1,KC
         QWSEDA(L,K)=0.0
      ENDDO
CGR      
        DO NL=1,NLWSEDA
          I=ICWSEDA(NL)
          J=JCWSEDA(NL)
          L=LIJ(I,J)
          TFACTOR=1.0
!          IF(IGRIDV.EQ.1) TFACTOR=GVCSCLP(L) !SCJ not GVC here
           NS=NWSESERA(NL)
           WSEMODEL=HP(L)+BELV(L)
           WSEOBSER=GI*PSERT(NS)
            IF(WSEOBSER.GT.RVAL)THEN
             QSSWSE=DELTI*TSWSEDA(NL)*DXYP(L)*(WSEOBSER-WSEMODEL)
             DO K=1,KC  !ADDED LOGIC TO ACCOUNT FOR DUPLICATE LOCATIONS greg rocheleau
              IF(QWSEDA(L,K).EQ.0) THEN
                QWSEDA(L,K)=DZC(K)*TFACTOR*QSSWSE
              ELSE
                IF (I.EQ.J)THEN
                   QWSEDA(L,K)=(QWSEDA(L,K)+DZC(K)*TFACTOR*QSSWSE)
     &              /REAL(2)
                ENDIF   
                IF(I.GT.J)THEN
                    QWSEDA(L,K)=QWSEDA(L,K)
                ENDIF    
                IF(J.GT.I)THEN
                    QWSEDA(L,K)=DZC(K)*TFACTOR*QSSWSE
                ENDIF
              ENDIF
            ENDDO
            QWSEASM(NDAYA,NL)=QWSEASM(NDAYA,NL)+QSSWSE
CGR          ELSE
CGR            DO K=1,KC
CGR              QWSEDA(L,K)=0.0
CGR            ENDDO
          ENDIF
C
        ENDDO
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C
C **  HORIZONTAL VELOCITY DATA ASSIMILATION
C
C     CALCULATIONS VECTOR BODY FORCES (HAVING KINEMATIC STRESS UNITS)
C     NECESSARY TO FORCE THE THE COMPUTED HORIZONTAL VELOCITY COMPONENTS
C     TOWARD THE OBSERVED HORIZONTAL VELOCITY COMPONENTS.  THE
C     BODY FORCES ARE INSERTED INTO THE EXPLICIT MOMENTUM EQUATION
C     TERMS (CALEXP, CALEXP2T, CALEXPGVC) ON THE NEXT TIME STEP
C     THE COMPUTED BODY FORCES (STRESSES) ARE WRITTEN TO THE EXTERNAL 
C     MODE VELOCITY TIME SERIES OUTPUT FILES
C
C----------------------------------------------------------------------C
C
C ASSIMILATE User Defined LAYERS
CGR***3/25/10
      IF(ISUVDA.LT.3)THEN
C*************************
CGR 4/14/10  Use adaptive nudging.  Weak during time of outflow.
      DO NL=1,NLUVDA
        TSUTMP(NL)=TSUUDA(NL)
        TSVTMP(NL)=TSVVDA(NL)
      ENDDO
      DO NL=1,NLUVDA
        I=ICUVDA(NL)
        J=JCUVDA(NL)
        L=LIJ(I,J)
CGR NOTE: TRY REDUCING RELAXATION IF (VHDX(LN,KC)/VHDX(L,KC)).LT.0 (I.E. DIFFERENT DIRECTIONS)        
CGR COMPUTE NORMAL PROPAGATION        
      IF (NORMDIR(NL).GT.0.OR.NORMDIR(NL).LT.0)THEN
CGR NORMDIR SHOULD BE NEGATIVE IF NORMAL OPPOSITE DIRECTION OF VNORTH (I.E. SOUTHERN BOUNDARY)      
       IF (NORMDIR(NL).LT.0)THEN
             LN=LNC(L)
C            PHASESPD=-(VHDX1E(L)-VHDX2E(L))*DELTI*DYP(L)/(VHDX1E(LN)-VHDX1E(L))
C_v1
C            PHASESPD=VHDXE(LN)*NORMDIR(NL)
C_v5
            PHASESPD=VHDXE(LN)*NORMDIR(NL)/abs(NORMDIR(NL))
C            PHASESPD=VHDX(LN,KC)/VHDX(L,KC); 
        ENDIF    
        IF(NORMDIR(NL).GT.0)THEN
            LN=LSC(L)
C            PHASESPD=-(VHDX1E(L)-VHDX2E(L))*DELTI*DYP(L)/(VHDX1E(LN)-VHDX1E(L))
C_v1
C            PHASESPD=VHDXE(LN)*NORMDIR(NL)
C_v5
            PHASESPD=VHDXE(LN)*NORMDIR(NL)/ABS(NORMDIR(NL)) 
C              PHASESPD=VHDX(LN,KC)/VHDX(L,KC); 
        ENDIF     
CGR         
       IF(PHASESPD.LT.0)THEN
            TSUTMP(NL)=MIN(TSUTMP(NL),0.05)
            TSVTMP(NL)=MIN(TSVTMP(NL),0.05)
        ENDIF
       ENDIF 
      ENDDO
C-------CONSTRAIN ASSIMILATION TO DEPTHS ABOVE KCUVDA (CURRENTLY PLACED IN ainit.FOR)     
c      DO K=1,KC
c        TUVKC(K)=1.0
c      ENDDO
c      IF(KCUVDA.GT.0)THEN
c       DO K=1,KCUVDA
c        TUVKC(K)=(NINT(K)/NINT(KCUVDA))**6.0
c       ENDDO  
c      ENDIF
C************************       
      IF(ICALL_.EQ.1)THEN
      IF(ISUVDA.GE.1)THEN
C
      DO NL=1,NLUVDA
        I=ICUVDA(NL)
        J=JCUVDA(NL)
        L=LIJ(I,J)
        NS=NUVSERA(NL)
	  DO K=1,KC
	    IF(USERT(K,NS).GT.RVAL) UROTTMP(K)=USERT(K,NS)
	    IF(VSERT(K,NS).GT.RVAL) VROTTMP(K)=VSERT(K,NS)
        ENDDO
	  DO K=1,KC
	    IF(USERT(K,NS).GT.RVAL) 
     &      USERT(K,NS)=WINDSXX(L)*UROTTMP(K)+WINDSXY(L)*VROTTMP(K)
	    IF(VSERT(K,NS).GT.RVAL) 
     &      VSERT(K,NS)=WINDSYX(L)*UROTTMP(K)+WINDSYY(L)*VROTTMP(K)
        ENDDO
        IF(SWB(L).LT.0.5)THEN
	    DO K=1,KC
            USERT(K,NS)=0.
            VSERT(K,NS)=0.
          ENDDO
        ENDIF
      ENDDO
C********************************
CGR 4/16/10  ASSIMILATE ALL LAYERS    
      DO NL=1,NLUVDA 
C
        TMPIMP=1.-FSUVDA(NL)
        I=ICUVDA(NL)
        J=JCUVDA(NL)
        L=LIJ(I,J)
	  LN=LNC(L)
	  LE=L+1
        NS=NUVSERA(NL)
        IF(ISTL_.EQ.2)THEN
          DO K=1,KC
            IF(USERT(K,NS).GT.RVAL)THEN
              RFBODYFXI(L)=RFBODYFXI(L)+1.
              RFBODYFXI(LE)=RFBODYFXI(LE)+1.
 	        FBODYFX(L,K)=DELTI*TSUTMP(NL)*TUVKC(K)*HU(L)*(USERT(K,NS)
     &                                          -TMPIMP*U(L,K))
     &                    +FBODYFX(L,K)
              FBODYFXI(L)=FSUVDA(NL)*TSUTMP(NL)*TUVKC(K)
     &                   +FBODYFXI(L)
 	        FBODYFX(LE,K)=DELTI*TSUTMP(NL)*TUVKC(K)*HU(LE)*(USERT(K,NS)
     &                                          -TMPIMP*U(LE,K))
     &                     +FBODYFX(LE,K)
              FBODYFXI(LE)=FSUVDA(NL)*TSUTMP(NL)*TUVKC(K)
     &                    +FBODYFXI(LE)
            ELSE
	        FBODYFX(L,K)=0.0
              FBODYFXI(L)=0.0
	        FBODYFX(LE,K)=0.0
              FBODYFXI(LE)=0.0
            ENDIF
            IF(VSERT(K,NS).GT.RVAL)THEN
              RFBODYFYI(L)=RFBODYFYI(L)+1.
              RFBODYFYI(LN)=RFBODYFYI(LN)+1.
 	        FBODYFY(L,K)=DELTI*TSVTMP(NL)*TUVKC(K)*HV(L)*(VSERT(K,NS)
     &                                          -TMPIMP*V(L,K))
     &                    +FBODYFY(L,K)
              FBODYFYI(L)=FSUVDA(NL)*TSVTMP(NL)*TUVKC(K)
     &                   +FBODYFYI(L)
 	        FBODYFY(LN,K)=DELTI*TSVTMP(NL)*TUVKC(K)*HV(LN)*(VSERT(K,NS)
     &                                          -TMPIMP*V(LN,K))
     &                     +FBODYFY(LN,K)
              FBODYFYI(LN)=FSUVDA(NL)*TSVTMP(NL)*TUVKC(K)
     &                    +FBODYFYI(LN)
            ELSE
	        FBODYFY(L,K)=0.0
              FBODYFYI(L)=0.0
	        FBODYFY(LN,K)=0.0
              FBODYFYI(LN)=0.0
            ENDIF
          ENDDO
        ELSE
          DO K=1,KC
            IF(USERT(K,NS).GT.RVAL)THEN
              RFBODYFXI(L)=RFBODYFXI(L)+1.
              RFBODYFXI(LE)=RFBODYFXI(LE)+1.
 	        FBODYFX(L,K)=DELTI*TSUTMP(NL)*TUVKC(K)*HU(L)*(USERT(K,NS)
     &                                          -TMPIMP*U1(L,K))
     &                    +FBODYFX(L,K)
              FBODYFXI(L)=FSUVDA(NL)*TSUTMP(NL)*TUVKC(K)
     &                   +FBODYFXI(L)
 	        FBODYFX(LE,K)=DELTI*TSUTMP(NL)*TUVKC(K)*HU(LE)*(USERT(K,NS)
     &                                          -TMPIMP*U1(LE,K))
     &                     +FBODYFX(LE,K)
              FBODYFXI(LE)=FSUVDA(NL)*TSUTMP(NL)*TUVKC(K)
     &                    +FBODYFXI(LE)
            ELSE
	        FBODYFX(L,K)=0.0
              FBODYFXI(L)=0.0
	        FBODYFX(LE,K)=0.0
              FBODYFXI(LE)=0.0
            ENDIF
            IF(VSERT(K,NS).GT.RVAL)THEN
              RFBODYFYI(L)=RFBODYFYI(L)+1.
              RFBODYFYI(LN)=RFBODYFYI(LN)+1.
 	        FBODYFY(L,K)=DELTI*TSVTMP(NL)*TUVKC(K)*HV(L)*(VSERT(K,NS)
     &                                          -TMPIMP*V1(L,K))
     &                    +FBODYFY(L,K)
              FBODYFYI(L)=FSUVDA(NL)*TSVTMP(NL)*TUVKC(K)
     &                   +FBODYFYI(L)
 	        FBODYFY(LN,K)=DELTI*TSVTMP(NL)*TUVKC(K)*HV(LN)*(VSERT(K,NS)
     &                                          -TMPIMP*V1(LN,K))
     &                     +FBODYFY(LN,K)
              FBODYFYI(LN)=FSUVDA(NL)*TSVTMP(NL)*TUVKC(K)
     &                    +FBODYFYI(LN)
            ELSE
	        FBODYFY(L,K)=0.0
              FBODYFYI(L)=0.0
	        FBODYFY(LN,K)=0.0
              FBODYFYI(LN)=0.0
            ENDIF
          ENDDO
CGR          
        ENDIF       
       ENDDO
C****************************************************
      DO K=1,KC
      DO L=1,LC
	  IF(RFBODYFXI(L).GT.0.0) FBODYFX(L,K)=FBODYFX(L,K)/RFBODYFXI(L)
        IF(RFBODYFYI(L).GT.0.0) FBODYFY(L,K)=FBODYFY(L,K)/RFBODYFYI(L)
      ENDDO
      ENDDO
C
      DO L=1,LC
	  IF(RFBODYFXI(L).GT.0.0) FBODYFXI(L)=FBODYFXI(L)/RFBODYFXI(L)
        IF(RFBODYFYI(L).GT.0.0) FBODYFYI(L)=FBODYFYI(L)/RFBODYFYI(L)
      ENDDO
C     
      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** CONSTRAIN ASSIMILATION TO DEPTH AVERAGE
C
      IF(ICALL_.EQ.1)THEN
      IF(ISUVDA.EQ.1)THEN
C
      IF(IGRIDV.EQ.0)THEN
        DO L=1,LC
          AVGX=0.0
          AVGY=0.0
          DO K=1,KC
            AVGX=AVGX+DZC(K)*FBODYFX(L,K)
	      AVGY=AVGY+DZC(K)*FBODYFY(L,K)
          ENDDO
          DO K=1,KC
            FBODYFX(L,K)=AVGX
	      FBODYFY(L,K)=AVGY
          ENDDO
        ENDDO
      ENDIF
C
      !IF(IGRIDV.EQ.1)THEN !SCJ removed because we are not using GVC
      !  DO L=1,LC
      !    AVGX=0.0
      !    AVGY=0.0
      !    DO K=1,KC
      !      AVGX=AVGX+DZC(K)*GVCSCLU(L)*FBODYFX(L,K)
	     ! AVGY=AVGY+DZC(K)*GVCSCLV(L)*FBODYFY(L,K)
      !    ENDDO
      !    DO K=KGVCU(L),KC
      !      FBODYFX(L,K)=AVGX
      !    ENDDO
      !    DO K=KGVCV(L),KC
	     ! FBODYFY(L,K)=AVGY
      !    ENDDO
      !  ENDDO
      !ENDIF
C
      ENDIF
	ENDIF
C
C**********************************************************************C
C
C **  HORIZONTAL VELOCITY DATA ASSIMILATION
C
C      DIAGNOSTIC FOR OUTPUT
C
C----------------------------------------------------------------------C
C
C ** ASSIMILATE ALL LAYERS
C
      IF(ICALL_.EQ.2)THEN
      IF(ISUVDA.GE.1)THEN
C
      DO NL=1,NLUVDA
C
        TMPIMP=1.-FSUVDA(NL)
        I=ICUVDA(NL)
        J=JCUVDA(NL)
        L=LIJ(I,J)
	  LN=LNC(L)
	  LE=L+1
        NS=NUVSERA(NL)
        IF(ISTL_.EQ.2)THEN
          DO K=1,KC
            IF(USERT(K,NS).GT.RVAL)THEN
 	        FBODYFX(L,K)=DELTI*TSUTMP(NL)*TUVKC(K)*HU(L)*(USERT(K,NS)
     &                       -FSUVDA(NL)*U(L,K)-TMPIMP*U1(L,K))
            ELSE
	        FBODYFX(L,K)=0.0
            ENDIF
            IF(VSERT(K,NS).GT.RVAL)THEN
 	        FBODYFY(L,K)=DELTI*TSVTMP(NL)*TUVKC(K)*HV(L)*(VSERT(K,NS)
     &                       -FSUVDA(NL)*V(L,K)-TMPIMP*V1(L,K))
            ELSE
	        FBODYFY(L,K)=0.0
            ENDIF
          ENDDO
        ELSE
          DO K=1,KC
            IF(USERT(K,NS).GT.RVAL)THEN
 	        FBODYFX(L,K)=DELTI*TSUTMP(NL)*TUVKC(K)*HU(L)*(USERT(K,NS)
     &                       -FSUVDA(NL)*U(L,K)-TMPIMP*U2(L,K))
            ELSE
	        FBODYFX(L,K)=0.0
            ENDIF
            IF(VSERT(K,NS).GT.RVAL)THEN
 	        FBODYFY(L,K)=DELTI*TSVTMP(NL)*TUVKC(K)*HV(L)*(VSERT(K,NS)
     &                       -FSUVDA(NL)*V(L,K)-TMPIMP*V2(L,K))
            ELSE
	        FBODYFY(L,K)=0.0
            ENDIF
          ENDDO
C
       ENDIF
       ENDDO
C
      ENDIF
      ENDIF
C
C----------------------------------------------------------------------C
C
C ** CONSTRAIN ASSIMILATION TO DEPTH AVERAGE
C
      IF(ICALL_.EQ.2)THEN
      IF(ISUVDA.EQ.1)THEN
C
      IF(IGRIDV.EQ.0)THEN
        DO NL=1,NLUVDA
          I=ICUVDA(NL)
          J=JCUVDA(NL)
          L=LIJ(I,J)
	    LN=LNC(L)
	    LE=L+1
          AVGX=0.0
          AVGY=0.0
          DO K=1,KC
            AVGX=AVGX+DZC(K)*FBODYFX(L,K)
	      AVGY=AVGY+DZC(K)*FBODYFY(L,K)
          ENDDO
          DO K=1,KC
            FBODYFX(L,K)=AVGX
	      FBODYFY(L,K)=AVGY
          ENDDO
        ENDDO
      ENDIF
C
!      IF(IGRIDV.EQ.1)THEN
!        DO NL=1,NLUVDA
!          I=ICUVDA(NL)
!          J=JCUVDA(NL)
!          L=LIJ(I,J)
!	    LN=LNC(L)
!	    LE=L+1
!          AVGX=0.0
!          AVGY=0.0
!          DO K=1,KC
!            AVGX=AVGX+DZC(K)*GVCSCLU(L)*FBODYFX(L,K)
!	      AVGY=AVGY+DZC(K)*GVCSCLV(L)*FBODYFY(L,K)
!          ENDDO
!          DO K=KGVCU(L),KC
!            FBODYFX(L,K)=AVGX
!          ENDDO
!          DO K=KGVCV(L),KC
!	      FBODYFY(L,K)=AVGY
!          ENDDO
!        ENDDO
!      ENDIF
C
      ENDIF
	ENDIF
C
C**********************************************************************C
CGR****3/25/10
      ENDIF
C      
      RETURN
      END
