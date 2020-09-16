      SUBROUTINE WQ3D(ISTL_,IS2TL_)  
C  
C  CONTROL SUBROUTINE FOR WATER QUALITY MODEL  
C  ORGINALLY CODED BY K.-Y. PARK  
C  OPTIMIZED AND MODIFIED BY J. M. HAMRICK  
C CHANGE RECORD  
C  
C     Merged SNL and DS-INTL
      USE GLOBAL  
#ifdef key_ncdf
      USE IOM
#endif
      REAL,SAVE :: DAYNEXT
      REAL,SAVE :: SUNDAY1, SUNDAY2, SUNSOL1, SUNSOL2
      REAL,SAVE :: SUNFRC1, SUNFRC2
      INTEGER*4,SAVE :: M
      
      DATA IWQTAGR,IWQTSTL/2*0/  
      DATA ISMTICI/0/  
      IF(ETIMEDAY.LE.(DTWQ+1.E-8))THEN
        DAYNEXT=FLOAT(INT(TIMEDAY))+1.
      ENDIF

      ! *** PMC - NEW IMPLEMENTATION TO USE DAILY (FROM HOURLY) SOLAR RADIATION FOR ALGAL GROWTH
      IF(ITNWQ.EQ.0.AND.IWQSUN.GT.1.AND.NASER.GT.0)THEN
        ! *** BUILD THE DAILY AVERAGE SOLAR RADIATION FROM THE ASER DATA
        SUNDAY1 = TIMEDAY+0.5
        SUNDAY2 = DAYNEXT+0.5
        
        ! *** FIND 1ST POINT
        M = 1
        DO WHILE (TASER(M,1).LT.SUNDAY1-0.5)
          M = M+1
        END DO
        
        ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
        M1 = 0
        M2 = 0
        SUNSOL1 = 0.0
        DO WHILE (TASER(M,1).LT.SUNDAY1+0.5)
          M1 = M1+1
          IF(SOLSWR(M,1).GT.0.)THEN
            M2 = M2+1
            SUNSOL1=SUNSOL1+SOLSWR(M,1)
          ENDIF
          M = M+1
        END DO
        IF(M1.GT.0)THEN
          SUNFRC1=FLOAT(M2)/FLOAT(M1)
          SUNSOL1=SUNSOL1/FLOAT(M1)
        ELSE
          SUNFRC1=1.0
        ENDIF
        
        ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
        M1 = 0
        M2 = 0
        SUNSOL2 = 0.
        DO WHILE (TASER(M,1).LT.SUNDAY2+0.5)
          M1 = M1+1
          IF(SOLSWR(M,1).GT.0.)THEN
            M2 = M2+1
            SUNSOL2=SUNSOL2+SOLSWR(M,1)
          ENDIF
          M = M+1
        END DO
        IF(M1.GT.0)THEN
          SUNFRC2=FLOAT(M2)/FLOAT(M1)
          SUNSOL2=SUNSOL2/FLOAT(M1)
        ELSE
          SUNFRC2=1.
        ENDIF
      ENDIF
C  
C **  READ INITIAL CONDITIONS  
C  
      IF(IWQICI.EQ.1) CALL RWQICI 
C  
C **  READ TIME/SPACE VARYING ALGAE PARAMETERS  
C  
      IF(IWQAGR.EQ.1 .AND. ITNWQ.EQ.IWQTAGR) CALL RWQAGR(IWQTAGR)  
C  
C **  READ TIME/SPACE VARYING SETTLING VELOCITIES  
C  
      IF(IWQSTL.EQ.1 .AND. ITNWQ.EQ.IWQTSTL) CALL RWQSTL(IWQTSTL)  
C  
C *** READ BENTHIC FLUX IF REQUIRED  
C *** CALL SPATIALLY AND TIME VARYING BENTHIC FLUX HERE.  ONLY CALL RWQBEN2  
C *** IF SIMULATION TIME IS >= THE NEXT TIME IN THE BENTHIC FILE.  
C  
      IF(IWQBEN .EQ. 2)THEN  
        IF(ISDYNSTP.EQ.0)THEN  
          TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.  
        ELSE  
          TIMTMP=TIMESEC/86400.  
        ENDIF  
        IF(TIMTMP .GE. BENDAY)THEN  
          CALL RWQBEN2(TIMTMP)  
        ENDIF  
      ENDIF  
C  
C **  UPDATE POINT SOURCE LOADINGS  
C  
      IF(IWQPSL.EQ.1)THEN
        CALL RWQPSL  
      ELSEIF(IWQPSL.GE.2) THEN !GR modified to be G.E.2 for option IWQPSL.EQ.3
        CALL CALCSER(ISTL_)
      ENDIF
C
      CALL RWQATM  
C  
C **  READ SEDIMENT MODEL INITIAL CONDITION  
C  
      IF(IWQBEN.EQ.1)THEN  
        IF(ISMICI.EQ.1 .AND. ITNWQ.EQ.ISMTICI) CALL RSMICI(ISMTICI)  
      ENDIF  
C  
C **  UPDATE OLD CONCENTRATIONS  
C   FOLLOWING THE CALL TO CALWQC MINUS OLD D.O. BEFORE THE CALL).  
C   FIRST SUBTRACT THE OLD D.O. HERE:  
C  
      IF(ISMTSB.LT.ISMTSE)THEN  
        DO K=1,KC  
          DO L=2,LA  
            XMRM = WQV(L,K,19)*DTWQ*DZC(K)*HP(L)  
            XDOTRN(L,K) = XDOTRN(L,K) - XMRM  
            XDOALL(L,K) = XDOALL(L,K) - XMRM  
          ENDDO  
        ENDDO  
      ENDIF
C  
C **  CALCULATE PHYSICAL TRANSPORT  
C **  WQV(L,K,NW) SENT TO PHYSICAL TRANSPORT AND TRANSPORTED  
C **  VALUE RETURNED IN WQV(L,K,NW)  
C  
      CALL CALWQC(ISTL_,IS2TL_) !transports (advects/disperses) WQV
C  
C   FOLLOWING THE CALL TO CALWQC MINUS OLD D.O. BEFORE THE CALL).  
C   NOW ADD THE NEW D.O. HERE:  
C  
      IF(ISMTSB.LT.ISMTSE)THEN  
        DO K=1,KC  
          DO L=2,LA  
            XMRM = WQV(L,K,19)*DTWQ*DZC(K)*HP(L)  
            XDOTRN(L,K) = XDOTRN(L,K) + XMRM  
            XDOALL(L,K) = XDOALL(L,K) + XMRM  
          ENDDO  
        ENDDO  
      ENDIF
C  
C **  UPDATE WATER COLUMN KINETICS AND SEDIMENT MODEL  
C **  OVER LONGER TIME INTERVALS THAN PHYSICAL TRANSPORT  
C **  IF NWQKDPT .GT. 1  
C  
      NWQKCNT=NWQKCNT+1  
      IF(ITNWQ.EQ.0.OR.NWQKCNT.EQ.NWQKDPT)THEN  
        !IF(ITNWQ.NE.0)NWQKCNT=0   PMC
        NWQKCNT=0
        ! **  UPDATE SOLAR RADIATION INTENSITY  
        !   WQI1 = SOLAR RADIATION ON PREVIOUS DAY  
        !   WQI2 = SOLAR RADIATION TWO DAYS AGO  
        !   WQI3 = SOLAR RADIATION THREE DAYS AGO  
        ! ***  UPDATE OCCURS ONLY WHEN THE SIMULATION DAY CHANGES.  
        IF(TIMEDAY.GT.DAYNEXT)THEN  ! *** DSLLC: FORCE A SOLAR DAY UPDATE
          WQI3 = WQI2  
          WQI2 = WQI1  
          WQI1 = WQI0OPT  
          IF(IWQSUN.GT.0)WQI0OPT = 0.0  
          DAYNEXT=DAYNEXT+1.
        ENDIF
        
        IF(IWQSUN.GT.1)THEN  
          IF(TIMEDAY.GT.SUNDAY2)THEN
            ! *** BUILD THE DAILY AVERAGE SOLAR RADIATION FROM THE ASER DATA
            SUNDAY1 = SUNDAY2
            SUNSOL1 = SUNSOL2
            SUNFRC1 = SUNFRC2
          ! *** FIND 1ST POINT
            !IF(M==0)M = 1 !If NASER = 0, the first time through M=0 so TASER(0,1) is not defined (this means there is an error in the INP files)
            ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION       
            M1 = 0
            M2 = 0
            SUNSOL2 = 0.
            SUNDAY2 = SUNDAY2+1.
            DO WHILE (TASER(M,1).LT.SUNDAY2+0.5)
              M1 = M1+1
              IF(SOLSWR(M,1).GT.0.)THEN
                M2 = M2+1
                SUNSOL2=SUNSOL2+SOLSWR(M,1)
              ENDIF
              M = M+1
            END DO
            IF(M1.GT.0)THEN
              SUNFRC2=FLOAT(M2)/FLOAT(M1)
              SUNSOL2=SUNSOL2/FLOAT(M1)
            ELSE
              SUNFRC2=1.
            ENDIF
            
          ENDIF
        ENDIF  
  
        ! **  READ SOLAR RADIATION INTENSITY AND DAYLIGHT LENGTH  
        ! NOTE: IWQSUN=1 CALLS SUBROUTINE RWQSUN WHICH READS THE DAILY  
        !                SOLAR RADIATION DATA FROM FILE SUNDAY.INP WHICH  
        !                ARE IN UNITS OF LANGLEYS/DAY.  
        !       IWQSUN=2 USES THE HOURLY SOLAR RADIATION DATA FROM ASER.INP  
        !                COUPLED WITH THE COMPUTED OPTIMAL DAILY LIGHT TO
        !                LIMIT ALGAL GROWTH.
        !       IWQSUN=3 USES THE DAILY AVERAGE SOLAR RADIATION DATA COMPUTED 
        !                FROM THE HOURLY ASER.INP AND THE COMPUTED OPTIMAL DAILY
        !                LIGHT TO LIMIT ALGAL GROWTH.
        !       IWQSUN>1 USES THE DAILY AVERAGE SOLAR RADIATION DATA COMPUTED 
        !                FROM THE HOURLY ASER.INP DATA.  CONVERTS WATTS/M**2 TO
        !                LANGLEYS/DAY USING 2.065.  COMPUTES THE FRACTION OF
        !                DAYLIGHT AND ADJUSTS FOR PHOTOSYNTHETIC ACTIVE RADIATION BY 
        !                PARADJ (~0.43) 
        !  
        IF(IWQSUN.EQ.0)THEN
          WQI0OPT = WQI0
        ELSEIF(IWQSUN.EQ.1)THEN  
          CALL RWQSUN  
          WQI0=SOLSRDT  
          WQFD=SOLFRDT  
          ! *** OPTIMAL SOLAR RADIATION IS ALWAYS UPDATED BASED ON DAY AVERAGED
          WQI0OPT = MAX(WQI0OPT, WQI0)
        ELSEIF(IWQSUN.GT.1)THEN
          RATIO = (TIMEDAY-SUNDAY1)
          SOLARAVG = RATIO*(SUNSOL2-SUNSOL1)+SUNSOL1
          WQFD=RATIO*(SUNFRC2-SUNFRC1)+SUNFRC1

          ! *** SOLAR RADIATION IN LANGLEYS/DAY
          WQI0 = PARADJ*2.065*SOLARAVG  

          IF(IWQSUN.EQ.2)THEN
            ! *** OPTIMAL SOLAR RADIATION IS ALWAYS UPDATED BASED ON DAY AVERAGED
            WQI0OPT = MAX(WQI0OPT, WQI0*0.85)

            IF(NASER.GT.1.OR.USESHADE)THEN  
              SOLARAVG=0.  
              DO L=2,LA  
                SOLARAVG=SOLARAVG+SOLSWRT(L)  
              ENDDO  
              SOLARAVG=SOLARAVG/FLOAT(LA-1)
            ELSE
              ! *** Spatially Constant Atmospheric Parameters
              SOLARAVG=SOLSWRT(2)
            ENDIF  
            ! *** SOLAR RADIATION IN LANGLEYS/DAY
            WQI0 = PARADJ*2.065*SOLARAVG  
            WQFD=1.  
          ELSE
            ! *** OPTIMAL SOLAR RADIATION IS ALWAYS UPDATED BASED ON DAY AVERAGED
            WQI0OPT = MAX(WQI0OPT, WQI0)  
          ENDIF
        ENDIF  
C  
C **  LOAD WQV INTO WQVO FOR REACTION CALCULATION  
C  
        NMALG=0  
        IF(IDNOTRVA.GT.0) NMALG=1  
        DO NW=1,NWQV+NMALG  
          IF(ISTRWQ(NW).NE.0)THEN  
!            DO K=1,KC  
!              DO L=2,LA  
!                WQVO(L,K,NW)=WQV(L,K,NW)  
!              ENDDO  
!            ENDDO
            WQVO(2:LA,1:KC,NW)=WQV(2:LA,1:KC,NW)
          ENDIF  
        ENDDO  
C  
C **    CALCULATE KINETIC SOURCES AND SINKS  
C  
        TTMP=SECNDS(SECNDS_ZERO)  
        IF(ISWQLVL.EQ.0) CALL WQSKE0  
        IF(ISWQLVL.EQ.1) CALL WQSKE1  
        IF(ISWQLVL.EQ.2) CALL WQSKE2  
        IF(ISWQLVL.EQ.3) CALL WQSKE3  
        IF(ISWQLVL.EQ.4) CALL WQSKE4  
        TWQKIN=TWQKIN+SECNDS(TTMP)  
C  
C **    DIAGNOSE NEGATIVE CONCENTRATIONS  
C  
        IF(IWQNC.EQ.1)CALL WWQNC  
C  
C **    WRITE TIME SERIES  
C
C *** Write of WQ spatial maps depends on netCDF library. If this is not
C *** available then it will not write spatial maps.  A user-defined function
C *** could be temporarily introduced that dumps WQ variables to ASCII file
C *** to investigate further, have a look at WQ_NC_WRITE in netcdf_iom.F90 file
        IF(ITNWQ.GT.0 .AND. MOD(ITNWQ,IWQTSDT).EQ.0 .AND. ISVPH.EQ.1)THEN
#ifdef key_ncdf
          CALL WQ_NC_WRITE
#endif
        ENDIF
        IF(ITNWQ.GE.IWQTSB .AND. ITNWQ.LE.IWQTSE.AND.IWQTSE.GT.0)THEN  
          IF(MOD(ITNWQ,IWQTSDT).EQ.0) CALL WWQTS
C  
          CALL WWQTSBIN  
C  
        ENDIF  
C  
C **    CALL SEDIMENT DIAGENSIS MODEL  
C  
        IF(IWQBEN.EQ.1)THEN  
          TTMP=SECNDS(SECNDS_ZERO)  
          CALL SMMBE  
          TWQSED=TWQSED+SECNDS(TTMP)  
          IF(ISMTS.GE.1)THEN  
C  
C **      WRITE SEDIMENT MODEL TIME SERIES  
C  
            IF(ITNWQ.GE.ISMTSB .AND. ITNWQ.LE.ISMTSE)THEN  
              IF(MOD(ITNWQ,ISMTSDT).EQ.0) CALL WSMTS
            ENDIF  
          ENDIF  
C  
C **      WRITE SEDIMENT MODEL FLUXES TO BINARY FILE:  
C  
          IF(ITNWQ.GE.ISMTSB .AND. ITNWQ.LE.ISMTSE)THEN  
            CALL WSMTSBIN  
          ENDIF  
        ENDIF  
      ENDIF  
C  
C **  UPDATE TIME IN DAYS  
C  
      ITNWQ = ITNWQ + 2  !Integer time number in days
C  
C **  ENDIF ON KINETIC AND SEDIMENT UPDATE  
C **  INSERT TIME CALL  
C **  WRITE RESTART FILES  
C  
      RETURN  
      END  

