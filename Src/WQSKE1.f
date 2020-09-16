      SUBROUTINE WQSKE1  
C  
C  ORGINALLY CODED BY K.-Y. PARK  
C  OPTIMIZED AND MODIFIED BY J.M. HAMRICK 
C 
C CHANGE RECORD  
C
C     MAJOR REWRITE BY PAUL M. CRAIG  JANUARY 12, 2006
      USE GLOBAL  

      IMPLICIT NONE
      
      INTEGER NQ,NS,IZ,IMWQZ,NSTPTMP,M
      INTEGER K,L,LE,LN
      
      REAL WQAVGIO,RMULTMP,TIME,RLIGHT1,RLIGHT2
      REAL WQGNC,WQGND,WQGNG,WQGNM,WQGPM,WQF1NM,WQGPC,WQGPD,WQGPG
      REAL WQF1NC,WQF1ND,WQF1NG,WQKESS,XMRM,YMRM,WQTT1
      
      REAL WQF2IC,WQF2ID,WQF2IG,SADWQ,WQGSD,WQF2IM
      REAL UMRM,VMRM,WQVEL,WQLVF,WQF4SC,WQKDOC,WQKHP,WQTTS
      REAL WQKHN,WQTTM,TVAL1,TVAL2,TVAL3,TVAL4,TVAL5
      REAL RLNSAT1,RLNSAT2,XNUMER,XDENOM,WQLDF,WQTTC,WQTTD,WQTTG
      REAL WINDREA,WQWREA,WQVREA,WQA1C,WQVA1C,WQR1C,WQA2D
      REAL WQR2D,WQA3G,WQR3G,WQB4,WQA4,WQR4,WQC5,WQA5,WQR5
      REAL WQD6,WQA6C,WQA6D,WQA6G,WQA6,WQA6M,WQR6
      REAL WQE7,WQA7C,WQA7D,WQA7G,WQA7,WQR7
      REAL WQF8,WQA8C,WQA8D,WQA8G,WQA8,WQR8
      REAL WQF9,WQA9C,WQA9D,WQA9G,WQA9,WQR9
      REAL WQA10C,WQA10D,WQA10G,WQKKL
      !REAL WQR10
      REAL WQI11,WQA11C,WQA11D,WQA11G,WQA11,WQR11
      REAL WQJ12,WQA12C,WQA12D,WQA12G,WQA12,WQR12
      REAL WQF13,WQA13C,WQA13D,WQA13G,WQA13,WQR13
      REAL WQR14,WQF14,WQA14C,WQA14D,WQA14G,WQA14
      REAL WQR15,WQA15C,WQA15D,WQA15G,WQA15,WQB15
      REAL WQM16,WQA16D,WQR16,WQR17,WQR18
      REAL TMP19,TEMFAC,DTWQxH,DTWQxH2,WQA19C,WQA19D,WQA19G
      REAL WQA19,WQA19A,WQSUM,WQRea,WQPOC,WQDOC,WQNH3,WQCOD
      REAL WQT20,WQR21,TIMTMP,WQTAMD
      REAL PPCDO,WQA22, WQA22C, WQA22D, WQA22G, WQCDDOC
      REAL WQCDREA, WQCDSUM
      REAL EXPA0,EXPA1 !VARIABLES FOR LIGHT EXTINCTION
      REAL WQGCO2M,WQGCO2C,WQGCO2G,WQGCO2D		!  CO2 Limitation Consts 
     
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DZCHP
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQISC
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQISD
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQISG
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQISM
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WQI0TOP
C
      ! ***  1) CHC - cyanobacteria
      ! ***  2) CHD - diatom algae
      ! ***  3) CHG - green algae
      ! ***  4) ROC - refractory particulate organic carbon
      ! ***  5) LOC - labile particulate organic carbon
      ! ***  6) DOC - dissolved organic carbon
      ! ***  7) ROP - refractory particulate organic phosphorus
      ! ***  8) LOP - labile particulate organic phosphorus
      ! ***  9) DOP - dissolved organic phosphorus
      ! *** 10) P4D - total phosphate
      ! *** 11) RON - refractory particulate organic nitrogen 23) macroalgae
      ! *** 12) LON - labile particulate organic nitrogen
      ! *** 13) DON - dissolved organic nitrogen
      ! *** 14) NHX - ammonia nitrogen
      ! *** 15) NOX - nitrate nitrogen
      ! *** 16) SUU - particulate biogenic silica
      ! *** 17) SAA - dissolved available silica
      ! *** 18) COD - chemical oxygen demand
      ! *** 19) DOX - dissolved oxygen
      ! *** 20) TAM - total active metal
      ! *** 21) FCB - fecal coliform bacteria
      ! *** 22) CO2 - dissolved carbon dioxide
      ! *** 23) macroalgae

      ! *** DTWQ - Water quality time step, which is typically in units of days
      ! *** DTWQO2 = DTWQ*0.5 

      ! *** WQCHL   = Chlorophyll a (ug/l)
	! *** WQCHLC  = carbon-to-chlorophyll ratio for cyanobacteria (mg C / ug Chl)
	! *** WQCHLD  = carbon-to-chlorophyll ratio for algae diatoms (mg C / ug Chl)
	! *** WQCHLG  = carbon-to-chlorophyll ratio for algae greens (mg C / ug Chl)
      ! *** WQKECHL = Light Extinction Coeff for CHLa (1/m per mg/l)
	! *** WQKETSS = Light Extinction Coeff for TSS (1/m per mg/l)
	! *** WQKEPOM = Light Extinction Coeff for POM (1/m per mg/l)
      ! *** WQKETOT(L,K) = Total Light Extinction
	! *** WQDOPG  = Optimal Depth for Growth - Green Algae

      ! *** RNH4WQ(L)  = Ammonia (for Current Layer)
      ! *** RNO3WQ(L)  = Nitrate (for Current Layer)
      ! *** PO4DWQ(L)  = Phosphate (for Current Layer)
      ! *** RNH4NO3(L) = Total Inorganic Nitrogen (for Current Layer)

      ! *** WQKHNG = Nitrogen half-saturation for Algae-Greens (mg/L)
      ! *** WQKHPG = Phosphorus half-saturation for Algae-Greens (mg/L)
	! *** WQKHCO2G = CO2 half-saturation for Algae-Greens (mg/L)

	! *** XLIMIG = Rate Limiting Factor - Light
	! *** XLIMTG = Rate Limiting Factor - Temperature   (Lookup Table: WQTDGG)
	! *** XLIMNG = Rate Limiting Factor - Nitrogen      (Local-WQGNG)
	! *** XLIMPG = Rate Limiting Factor - Phosphorus    (Local-WQGPG)
	! *** XLIMCO2G = Rate Limiting Factor - CO2    (Local-WQGPG) 
      ! *** WQF1NG = Rate Limiting Factor, Minimum of N & P

      ! *** WQPMG  = Maximum Growth Rate for Algae-Greens (1/d)
      ! *** WQPG   = Current Growth Rate for Algae-Greens (1/d)
      ! *** WQBMG  = Current Basal Metabolism Rate for Algae-Greens (1/d)
      ! *** WQPRG  = Current Predation Metabolism Rate for Algae-Greens (1/d)

      ! *** WQBMRG   = Basal Metabolism Rate for Algae-Greens (1/d)
      ! *** WQPRRG   = Predation Rate for Algae-Greens (1/d)
      ! *** WQTDRG   = Lookup Table for Temperature Rate Effect - Algae-Greens

      ! *** WQPC   = Final Net Growth Rate - Cyanobacteria
      ! *** WQPD   = Final Net Growth Rate - Diatoms Algae  
      ! *** WQPG   = Final Net Growth Rate - Green Algae
      ! *** WQPM   = Final Net Growth Rate - Macroalgae
      
      ! *** WQPNC = Preference for ammonium uptake - Cyanobacteria
      ! *** WQPND = Preference for ammonium uptake - Diatoms Algae
      ! *** WQPNG = Preference for ammonium uptake - Green Algae

      ! *** WQOBTOT  = Total Algal Biomass (mg/l)
      ! *** WQKRC    = Minimum Dissolution Rate of Refractory POC (1/day)
      ! *** WQKLC    = Minimum Dissolution Rate of Labile POC (1/day)
      ! *** WQKLCALG = Constant Refractory POC Dissolution Rate
      ! *** WQTDHDR  = Lookup Table for Temperature Rate Effect for Hydrolysis
      ! *** WQKRPC   = Current Dissolution Rate for POC

      ! *** WQI0   = SOLAR RADIATION for Current Time
      ! *** WQI1   = SOLAR RADIATION ON PREVIOUS DAY  
      ! *** WQI2   = SOLAR RADIATION TWO DAYS AGO  
      ! *** WQI3   = SOLAR RADIATION THREE DAYS AGO  

      ! *** WQKHR  = DOC Heterotrophic Respiration Rate

      ! *** WQWSSET = Water quality settling speed L;(L:1) is for top water layer; (L,2) is for lower water layers
      ! *** WQTTM   = Temporary concentration variable

      IF(.NOT.ALLOCATED(DZCHP))THEN
        ALLOCATE(DZCHP(LCM))   
        ALLOCATE(WQISC(LCM))  
        ALLOCATE(WQISD(LCM))  
        ALLOCATE(WQISG(LCM))  
        ALLOCATE(WQISM(LCM)) 
        ALLOCATE(WQI0TOP(LCM))  

        DZCHP=0.0
      ENDIF
C
      
      NS=1
      WQKESS=0.0
C    
C COMPUTE WQCHL,WQTAMP,WQPO4D,WQSAD AT A NEW TIME STEP: WQCHLX=1/WQCHLX  
C  
      ! *** Compute WQCHL (Chlorophyll) Using Algal Biomass & factors
      DO K=1,KC
        DO L=2,LA
          WQCHL(L,K)=WQV(L,K,1)*WQCHLC
     &              +WQV(L,K,2)*WQCHLD
     &              +WQV(L,K,3)*WQCHLG
        ENDDO
      ENDDO
C
C INITIALIZE SOLAR RADIATION AND OPTIMAL LIGHT
C
      ! *** INITIAL SOLAR RADIATION AT TOP OF SURFACE LAYER
      IF(USESHADE)THEN
        DO L=2,LA
          WQI0TOP(L)=WQI0 * PSHADE(L)
        ENDDO
      ELSE
        DO L=2,LA  
          WQI0TOP(L)=WQI0
        ENDDO  
      ENDIF
      ! ***  COMPUTE THE CURRENT OPTIMAL LIGHT INTENSITY
      IF(IWQSUN==2)THEN  
        WQAVGIO = WQCIA*WQI1 + WQCIB*WQI2 + WQCIC*WQI3  
      ELSE
        WQAVGIO = WQCIA*WQI0 + WQCIB*WQI1 + WQCIC*WQI2
      ENDIF 
      ! *** CORRECT TO AVERAGE SOLAR RADIATION DURING DAYLIGHT HOURS
 !     WQAVGIO = WQAVGIO / (WQFD + 1.E-18)											!!!!!!!!!!

      ! *** DZWQ=1/H (for a layer), VOLWQ=1/VOL (m^-3)
      DO K=KC,1,-1
        DO L=2,LA
          TWQ(L)=TEM(L,K)                 !layer temperature for WQ calcs
          SWQ(L)=MAX(SAL(L,K), 0.0)       !layer salinity for WQ calcs
          DZCHP(L)=DZC(K)*HP(L)           !layer thickness of a cell in meters
          DZWQ(L) = 1.0 / DZCHP(L)        !inverse layer thickness
          VOLWQ(L) = DZWQ(L) / DXYP(L) !inverse volume of each cell in a layer
          IMWQZT(L)=IWQZMAP(L,K)          !binary map for WQ calcs
        ENDDO    
        ! *** ZERO WQWPSL IF FLOWS ARE NEGATIVE.  THESE ARE HANDLED IN CALFQC (PMC)
        IF(IWQPSL/=2)THEN
          DO NQ=1,NQSIJ  
            IF((QSERCELL(K,NQ)+QSS(K,NQ))<=0.0)WQWPSL(LQS(NQ),K,1:NWQV)=0.0 ! *** ZERO THE FLUX
          ENDDO
        ENDIF
        IF(ISTRWQ(1)==1.OR.ISTRWQ(2)==1.OR.ISTRWQ(3)==1)THEN !EVALUATING THE RATE OF ALGAE LEAVING THE CELL THROUGH SETTLING OR FLOATING    
	    DO L=2,LA  
		  IF(WQWSC(IMWQZT(L))<0.0) THEN		!VB PERMITS CYANOBACTERIA TO FLOAT AND/OR SETTLE
		    IF(K==KC) THEN			
			  WQBCSET(L,1) = 0.0			!ALGAE AT THE WATER SURFACE CANT FLOAT INTO THE CELL ABOVE
              ELSE
		      WQBCSET(L,1) = -WQWSC(IMWQZT(L))*DZWQ(L)   ! *** CYANOBACTERIA	!NEEDS TO BE A POSITIVE QTY	
              ENDIF
		  ELSE
              WQBCSET(L,1) = WQWSC(IMWQZT(L))*DZWQ(L)   ! *** CYANOBACTERIA
		  ENDIF
		  IF(WQWSD(IMWQZT(L))<0.0) THEN		!VB PERMITS DIATOMS TO FLOAT AND/OR SETTLE
              IF(K==KC) THEN
			  WQBDSET(L,1) = 0.0
              ELSE
			  WQBDSET(L,1) = -WQWSD(IMWQZT(L))*DZWQ(L)   ! *** Diatoms	
              ENDIF
		  ELSE
		    WQBDSET(L,1) = WQWSD(IMWQZT(L))*DZWQ(L)   ! *** Diatoms
		  ENDIF
		  IF(WQWSG(IMWQZT(L))<0.0) THEN		!VB PERMITS GREEN ALGAE TO FLOAT AND/OR SETTLE
              IF(K==KC) THEN
			  WQBGSET(L,1) = 0.0
              ELSE
		  	  WQBGSET(L,1) = -WQWSG(IMWQZT(L))*DZWQ(L)   ! *** ALGAE	
              ENDIF
		  ELSE
		    WQBGSET(L,1) = WQWSG(IMWQZT(L))*DZWQ(L)   ! *** ALGAE
		  ENDIF
	    ENDDO
        ENDIF
        ! *** ZONE SPECIFIC SETTING VELOCITIES, (m/day)   
        IF(ISTRWQ(7)==1)WQRPSET(2:LA,1) = WQWSRP(IMWQZT(2:LA))*DZWQ(2:LA)  ! *** Refractory POM
        IF(ISTRWQ(8)==1)WQLPSET(2:LA,1) = WQWSLP(IMWQZT(2:LA))*DZWQ(2:LA)  ! *** Labile POM 
        ! *** SET SETTLING FOR TAM SORPTION: CURRENT LAYER  
        IF(IWQSRP==1)WQWSSET(2:LA,1) = WQWSS(IMWQZT(2:LA))*DZWQ(2:LA)  

        IF(K/=KC)IMWQZT1(2:LA)=IWQZMAP(2:LA,K+1)  
        IF(K/=1)IMWQZT2(2:LA)=IWQZMAP(2:LA,K-1)  
        
        IF(ISTRWQ(1)==1.OR.ISTRWQ(2)==1.OR.ISTRWQ(3)==1)THEN
          DO L=2,LA
		  IF(K/=KC)THEN 
	        IF(WQWSC(IMWQZT1(L))<0.0) THEN
			  WQBCSET(L,2) = 0.0  
              ELSE
			  WQBCSET(L,2) = WQWSC(IMWQZT1(L))*DZWQ(L)  
              ENDIF
              IF(WQWSD(IMWQZT1(L))<0.0) THEN
			  WQBDSET(L,2) = 0.0  
              ELSE
			  WQBDSET(L,2) = WQWSD(IMWQZT1(L))*DZWQ(L)  
              ENDIF
	  	    IF(WQWSG(IMWQZT1(L))<0.0) THEN
			  WQBGSET(L,2) = 0.0  
              ELSE
			  WQBGSET(L,2) = WQWSG(IMWQZT1(L))*DZWQ(L)  
              ENDIF
		  ENDIF
		  IF(K/=1)THEN
		    IF(WQWSC(IMWQZT2(L))<0.0)WQBCSET(L,2) = WQBCSET(L,2)-WQWSC(IMWQZT2(L))*DZWQ(L)  
              IF(WQWSD(IMWQZT2(L))<0.0)WQBDSET(L,2) = WQBDSET(L,2)-WQWSD(IMWQZT2(L))*DZWQ(L) 
              IF(WQWSG(IMWQZT2(L))<0.0)WQBGSET(L,2) = WQBGSET(L,2)-WQWSG(IMWQZT2(L))*DZWQ(L)  
		  ENDIF
          ENDDO
        ENDIF
        
        IF(K/=KC) THEN
          DO L=2,LA
            WQRPSET(L,2) = WQWSRP(IMWQZT1(L))*DZWQ(L)  
            WQLPSET(L,2) = WQWSLP(IMWQZT1(L))*DZWQ(L)
          ENDDO
        ENDIF
! *** SET SETTLING FOR TAM SORPTION: ONE LAYER UP
	  IF(IWQSRP==1.AND.K/=KC)THEN
          DO L=2,LA
            WQWSSET(L,2) = WQWSS(IMWQZT1(L))*DZWQ(L)
          ENDDO
        ENDIF  
C  
C FIND AN INDEX FOR LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY  
C  
		! *** DSLLC BEGIN BLOCK
        DO L=2,LA
          IWQT(L)=NINT((TWQ(L)-WQTDMIN)/WQTDINC)+1
        ENDDO  
        DO L=2,LA  
          IF(IWQT(L)<1 .OR. IWQT(L)>NWQTD)THEN  
            OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',POSITION='APPEND',STATUS='UNKNOWN')  
            WRITE(1,*)' *** ERROR IN WQSKE1:TEMPERATURE LOOKUP TABLE'
            WRITE(1,911) TIMEDAY, L, IL(L), JL(L), K, TWQ(L),TEM(L,K)  
            WRITE(6,600)IL(L),JL(L),K,TWQ(L),TEM(L,K)  
            IWQT(L)=MAX(IWQT(L),1)  
            IWQT(L)=MIN(IWQT(L),NWQTD)  
            CLOSE(1,STATUS='KEEP')
          ENDIF  
        ENDDO  
		! *** DSLLC END BLOCK
  600 FORMAT(' I,J,K,TWQ,TEM = ',3I5,2E13.4)  
  911 FORMAT('ERROR: TIME, L, I, J, K, TWQ(L),TEM(L,K) = ',  
     &    F10.5, I5, 2I4, I3, 2F10.4,/)  
  
        !C NOTE: MRM 04/29/99  ADDED ARRAYS TO KEEP TRACK OF  
        !C       NITROGEN, PHOSPHORUS, LIGHT, AND TEMPERATURE LIMITS  
        !C       FOR ALGAE GROWTH FOR CYANOBACTERIA, DIATOMS, GREENS,  
        !C       AND MACROALGAE.  THESE ARE THE ARRAYS:  
        !C        XLIMNX(L,K) = NITROGEN    LIMITATION FOR ALGAE GROUP X  
        !C        XLIMPX(L,K) = PHOSPHORUS  LIMITATION FOR ALGAE GROUP X  
        !C        XLIMIX(L,K) = LIGHT       LIMITATION FOR ALGAE GROUP X  
        !C        XLIMTX(L,K) = TEMPERATURE LIMITATION FOR ALGAE GROUP X  
       
        ! *** BEGIN HORIZONTAL LOOP FOR ALGAE PARMETERS
        DO L=2,LA
          RNH4WQ(L) = MAX (WQV(L,K,14), 0.0)     ! *** Ammonia
          RNO3WQ(L) = MAX (WQV(L,K,15), 0.0)     ! *** Nitrate
          PO4DWQ(L) = MAX (WQPO4D(L,K), 0.0)     ! *** Phosphate
          RNH4NO3(L) = RNH4WQ(L) + RNO3WQ(L)     ! *** Total Inorganic Nitrogen
          IF(ISTRWQ(22)==1)CO2WQ(L)  = MAX (WQV(L,K,22), 0.0)     ! *** CO2
        ENDDO
        DO L=2,LA  
          IF(ISTRWQ(1)==1.OR.ISTRWQ(2)==1.OR.ISTRWQ(3)==1)THEN !Microalgae?
            WQGNC = RNH4NO3(L) / (WQKHNC+RNH4NO3(L) + 1.E-18)  
            WQGND = RNH4NO3(L) / (WQKHND+RNH4NO3(L) + 1.E-18)  
            WQGNG = RNH4NO3(L) / (WQKHNG+RNH4NO3(L) + 1.E-18)
            WQGPC = PO4DWQ(L) / (WQKHPC+PO4DWQ(L) + 1.E-18)  
            WQGPD = PO4DWQ(L) / (WQKHPD+PO4DWQ(L) + 1.E-18)  
            WQGPG = PO4DWQ(L) / (WQKHPG+PO4DWQ(L) + 1.E-18)  
            WQGCO2C = CO2WQ(L) / (WQKHCO2C+CO2WQ(L) + 1.E-18)
            WQGCO2D = CO2WQ(L) / (WQKHCO2D+CO2WQ(L) + 1.E-18)
            WQGCO2G = CO2WQ(L) / (WQKHCO2G+CO2WQ(L) + 1.E-18)
            XLIMNC(L,K) = XLIMNC(L,K) + WQGNC  !Cyanobacteria nitrogen limitation
            XLIMND(L,K) = XLIMND(L,K) + WQGND  !Diatom        nitrogen limitation  
            XLIMNG(L,K) = XLIMNG(L,K) + WQGNG  !Green algae   nitrogen limitation
            XLIMPC(L,K) = XLIMPC(L,K) + WQGPC  !Cyanobacteria phosphorus limitation  
            XLIMPD(L,K) = XLIMPD(L,K) + WQGPD  !Diatom        phosphorus limitation   
            XLIMPG(L,K) = XLIMPG(L,K) + WQGPG  !Green algae   phosphorus limitation
            XLIMCO2C(L,K) = XLIMCO2C(L,K) + WQGCO2C  !Cyanobacteria CO2 limitation
            XLIMCO2D(L,K) = XLIMCO2D(L,K) + WQGCO2D  !Diatom        CO2 limitation  
            XLIMCO2G(L,K) = XLIMCO2G(L,K) + WQGCO2G  !Green algae   CO2 limitation
            XLIMTC(L,K) = XLIMTC(L,K) + WQTDGC(IWQT(L))  !Cyanobacteria temperature limitation  
            XLIMTD(L,K) = XLIMTD(L,K) + WQTDGD(IWQT(L))   !Diatom       temperature limitation 
            XLIMTG(L,K) = XLIMTG(L,K) + WQTDGG(IWQT(L))   !Green algae  temperature limitation 
            IF(ISTRWQ(22)>0)THEN 
              WQF1NC = MIN(WQGNC, WQGPC, WQGCO2C)			!Minimum of the N/P/CO2 Limit: Cyanobacteria		
	      ELSE
              WQF1NC = MIN(WQGNC, WQGPC)					!Minimum of the N/P     Limit: Cyanobacteria
            ENDIF
            IF(IWQSI==1)THEN !SILICA LIMITATION?
              SADWQ = MAX (WQSAD(L,K), 0.0)  
              WQGSD = SADWQ / (WQKHS+SADWQ+ 1.E-18)  
              IF(ISTRWQ(22)>0)THEN !CO2 LIMITATION?
                WQF1ND = MIN(WQGND, WQGPD, WQGSD, WQGCO2D) !Minimum of the N/P/Si/CO2 Limit: Diatoms
              ELSE
                WQF1ND = MIN(WQGND, WQGPD, WQGSD)          !Minimum of the N/P/Si     Limit: Diatoms
              ENDIF
            ELSEIF(IWQSI==0)THEN
              IF(ISTRWQ(22)>0)THEN !CO2 LIMITATON?    
                WQF1ND = MIN(WQGND, WQGPD, WQGCO2D)        !Minimum of the N/P/CO2    Limit: Diatoms
              ELSE
                WQF1ND = MIN(WQGND, WQGPD)                 !Minimum of the N/P        Limit: Diatom
              ENDIF
            ENDIF	
            IF(ISTRWQ(22)>0)THEN  
              WQF1NG = MIN(WQGNG, WQGPG, WQGCO2G)			!Minimum of the N/P/CO2   Limit: Greens			
	      ELSE
              WQF1NG = MIN(WQGNG, WQGPG)					!Minimum of the N/P       Limit: Greens	
            ENDIF
C  
C ALGAL BASAL METABOLISM & PREDATION  
C  
            WQBMC(L) = WQBMRC(IMWQZT(L)) * WQTDRC(IWQT(L))  
            WQPRC(L) = WQPRRC(IMWQZT(L)) * WQTDRC(IWQT(L))  
C  
C THE VARIABLE WQTDGP ADJUSTS PREDATION AND BASAL METABOLISM BASED ON A  
C LOWER/UPPER OPTIMUM TEMPERATURE FUNCTION.  THIS WILL ALLOW DIATOMS TO  
C BLOOM IN WINTER IF WQTDGP IS CLOSE TO ZERO.  
C  
            WQBMD(L)=WQBMRD(IMWQZT(L))*WQTDRD(IWQT(L))*WQTDGP(IWQT(L))  
            WQPRD(L)=WQPRRD(IMWQZT(L))*WQTDRD(IWQT(L))*WQTDGP(IWQT(L))  
            WQBMG(L) = WQBMRG(IMWQZT(L)) * WQTDRG(IWQT(L))  
            WQPRG(L) = WQPRRG(IMWQZT(L)) * WQTDRG(IWQT(L)) 
          ENDIF
!***For macroalgae defined in VEGE.INP IDNOTRVA>0
          IF(RMAC(L,K)>0.0)THEN !RMAC is the ratio of a layer occupied by macroalgae as calculated in CALTBXY
!***WQGNM is nitrate/ammonium limitation
            WQGNM = RNH4NO3(L) / (WQKHNM+RNH4NO3(L) + 1.E-18)
            MACLIM(L,K,4) = WQGNM !NO3/NH4 limitation saved
!***WQGPM is phosphate limitation
            WQGPM = PO4DWQ(L) / (WQKHPM+PO4DWQ(L) + 1.E-18)
            MACLIM(L,K,5) = WQGPM !PO4 limitation saved
!***WQGCO2M is CO2 limitation
!***WQF1NM is total nutrient limitation
            IF(ISTRWQ(22)>0)THEN								
              WQGCO2M = CO2WQ(L) / (WQKHCO2M+CO2WQ(L) + 1.E-18)	!CO2 macroalgae limitation
              MACLIM(L,K,6) = WQGCO2M !CO2 limitation saved
              WQF1NM = MIN(WQGNM, WQGPM, WQGCO2M)					!Minimum of the N/P/CO2 Limit: macroalgae
              XLIMCO2M(L,K) = XLIMCO2M(L,K) + WQGCO2M !CO2
            ELSE													
              WQF1NM = MIN(WQGNM, WQGPM)							!Minimum of the N/P     Limit: macroalgae
            ENDIF
            XLIMNM(L,K) = XLIMNM(L,K) + WQGNM  !Macroalgae nitrate/ammonium limitation
            XLIMPM(L,K) = XLIMPM(L,K) + WQGPM  !Macroalgae phosphorus limitation
            XLIMTM(L,K) = XLIMTM(L,K) + WQTDGM(IWQT(L))  !Macroalgae temperature limitation
            MACLIM(L,K,3) = WQTDGM(IWQT(L))    !Temperature limitation saved
!C BIOLOGICAL CARRYING CAPACITY LIMITATION
!C FIRST CONVERT FROM MACROALGAE FROM A CONCENTRATION (MG C/M3) TO A DENSITY (MG C/M2).
            IF(IWQVLIM>0)THEN !Is biomass carrying capacity a limitation?
              XMRM = WQV(L,K,IDNOTRVA)*DZCHP(L)  
              WQLDF = WQKBP(L) / (WQKBP(L) + XMRM) !Macroalgae biomass carrying capacity limitation
              XLIMDM(L,K) = XLIMDM(L,K) + WQLDF
            ELSE
              WQLDF=1.0 !No macroalgal biomass carrying capacity limitation on growth
            ENDIF
!C OPTION 1 FOR VELOCITY LIMITATION ASSUMES MACROALGAE GROWTH  
!C IS LIMITED AT LOW VELOCITIES DUE TO REDUCED AVAILABILITY OF  
!C NUTRIENTS REACHING THE ALGAE BIOMASS.  USES A MICHAELIS-MENTON  
!C TYPE OF EQUATION.  
            IF(IWQVLIM==1)THEN !Is macroalgal growth limited according to the Michaelis-Menton equation?
              LE=LEAST(L)
              LN=LNC(L)
              UMRM = 0.5*( U(L,K) + U(LE,K) )  
              VMRM = 0.5*( V(L,K) + V(LN,  K) )  
              WQVEL=SQRT(UMRM*UMRM + VMRM*VMRM)  
              IF(WQVEL>WQKMVMIN(L))THEN  
                WQLVF = WQVEL / (WQKMV(L) + WQVEL)  
              ELSE  
                WQLVF = WQKMVMIN(L) / (WQKMV(L) + WQKMVMIN(L))  
              ENDIF
!C OPTION 2 FOR VELOCITY LIMITATION APPLIES A FIVE-PARAMETER LOGISTIC  
!C FUNCTION THAT CAN BE ADJUSTED TO LIMIT MACROALGAE GROWTH FOR  
!C EITHER LOW OR HIGH (SCOUR) VELOCITIES.  IN STREAMS WITH LOW NUTRIENTS,  
!C THE LOW VELOCITY WILL LIKELY BE LIMITING SINCE AMPLE NUTRIENTS MAY  
!C NOT REACH THE ALGAE BIOMASS DUE TO REDUCED FLOW.  IN STREAMS WITH  
!C ABUNDANT NUTRIENTS, LOW VELOCITIES WILL NOT LIMIT MACROALGAE GROWTH,  
!C INSTEAD, HIGH VELOCITIES WILL LIKELY SCOUR THE MACROALGAE AND DETACH  
!C IT FROM THE SUBSTRATE.  
            ELSEIF(IWQVLIM==2)THEN !Is macroalgal growth limited according to a five-parameter logistic equation?
              XNUMER = WQKMVA(L) - WQKMVD(L)  
              XDENOM = 1.0 + (WQVEL/WQKMVC(L))**WQKMVB(L)  
              WQLVF = WQKMVD(L) + ( XNUMER / (XDENOM**WQKMVE(L)) )  
              XLIMVM(L,K) = XLIMVM(L,K) + WQLVF  !Macroalgae velocity limitation
            ELSE !No macroalgal velocity limitation on growth
              WQLVF=1.0
            ENDIF
! *** USE THE MORE SEVERELY LIMITING OF VELOCITY OR NUTRIENT FACTORS:  
            WQF1NM = MIN(WQLVF, WQF1NM)
! *** MACROALGAE BASAL METABOLISM AND PREDATION
            WQBMM(L) = WQBMRM(IMWQZT(L)) * WQTDRM(IWQT(L))
            WQPRM(L) = WQPRRM(IMWQZT(L)) * WQTDRM(IWQT(L))  
            MACLIM(L,K,7) = WQBMM(L) !Macroalgae basal metabolic rate saved
            MACLIM(L,K,8) = WQPRM(L) !Macroalgae predation rate saved
          ENDIF !End of macroalgae nutrient, biomass carrying capacity, and velocity limitation calculations
          ! *** IN C&C, F2IC=F2IC/FCYAN, FACTOR TO ALLOW CYANOBACTERIA MAT FORMATION
          ! *** COMPUTE TOTAL EXTINCTION COEFFICIENT
          ! *** LIGHT EXTINCTION (THIS WILL ALWAYS BE TRUE EXCEPT FOR IWQSUN=2)
          WQKESS=WQKEB(IMWQZT(L)) !Start with background light extinction
          IF(ISTRWQ(4)>0.OR.ISTRWQ(5)>0)WQKESS=WQKESS+WQKEPOM*(WQV(L,K,4)+WQV(L,K,5))*DZCHP(L) !Add any refractory and particulate organic carbon component
          IF(ISTRAN(6)>0.OR.ISTRAN(7)>0)WQKESS=WQKESS+WQKETSS*(SEDT(L,K)+SNDT(L,K))  *DZCHP(L) !Add any sediment component
          IF(ISTRWQ(1)==1.OR.ISTRWQ(2)==1.OR.ISTRWQ(3)==1)THEN !Add any Chlorophyll component
            IF(WQKECHL<0.0)THEN
              ! *** Compute Extinction Factor as a fn(Chla)
		    XMRM = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)
            ELSE
	        XMRM = WQKECHL*WQCHL(L,K)								
            ENDIF  
            WQKESS = WQKESS+XMRM
          ENDIF
          IF(RMAC(L,K)>0)WQKESS = WQKESS + WQKECHL*WQV(L,K,IDNOTRVA)*DZCHP(L) !Add any macroalgae component (may need its own light extinction, Ke, variable WQKEMAC) 
		IF(K==KC)THEN										!Specify surface solar radiation intensity
		  WQITOP(L,K) = WQI0TOP(L)
	    ELSE												!Calculate solar radiation intensity as a function of the layer above's solar radiation intensity
		  WQITOP(L,K) = WQITOP(L,K+1)*EXP(-WQKESS*DZCHP(L))
          ENDIF !SEE DiTORO ET AL. (1971, EQNS. (11)&(12))
! *** NOTE THAT LIGHT LIMITATION IS DUE TO EITHER Cyanobacteria, Diatoms, Green algae, or Macroalgae, BUT THESE ARE NOT ADDITIVE.
!     IN SYSTEMS WITH MORE THAT ONE MICROALGAE OR MICROALGAE PLUS MACROALGAE, THIS MUST BE REWRITTEN TO BE ADDITIVE.
          IF(WQI0>0.1.AND.(ISTRWQ(1)==1.OR.ISTRWQ(2)==1.OR.ISTRWQ(3)==1))THEN !If there is solar radiation and microalgae 
            ! *** OPTIMAL LIGHT INTENSITY AT OPTIMAL DEPTH
            IF(K==KC)THEN
              WQISC(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPC),WQISMIN)  
              WQISD(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPD),WQISMIN)  
              WQISG(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPG),WQISMIN)
            ENDIF
            ! *** LIGHT GROWTH-LIMITING FACTOR
		  EXPA0=EXP(-WQITOP(L,K)/WQISC(L))
		  EXPA1=EXP(-WQITOP(L,K)/WQISC(L)*EXP(-WQKESS*DZCHP(L)))
		  WQF2IC=EXP(1.0)*WQFD/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0) !Cyanobacteria light limitation
            XLIMIC(L,K) = XLIMIC(L,K) + WQF2IC
		  EXPA0=EXP(-WQITOP(L,K)/WQISD(L))
		  EXPA1=EXP(-WQITOP(L,K)/WQISD(L)*EXP(-WQKESS*DZCHP(L)))
		  WQF2ID=EXP(1.0)*WQFD/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0) !Diatom light limitation
            XLIMID(L,K) = XLIMID(L,K) + WQF2ID
		  EXPA0=EXP(-WQITOP(L,K)/WQISG(L))
		  EXPA1=EXP(-WQITOP(L,K)/WQISG(L)*EXP(-WQKESS*DZCHP(L)))
		  WQF2IG=EXP(1.0)*WQFD/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0) !Green algae light limitation
            XLIMIG(L,K) = XLIMIG(L,K) + WQF2IG
! *** Compute Microalgal Growth Rates due to Limitation Factors
            IF(IWQSTOX==1)THEN !Are cyanotoxins considered?
              WQF4SC = WQSTOX / (WQSTOX + SWQ(L)*SWQ(L)+1.E-12)  
              WQPC(L)=WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L))*WQF4SC  
            ELSE  
              WQPC(L) = WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L))  
            ENDIF  
            WQPD(L) = WQPMD(IMWQZT(L))*WQF1ND*WQF2ID*WQTDGD(IWQT(L))  
            WQPG(L) = WQPMG(IMWQZT(L))*WQF1NG*WQF2IG*WQTDGG(IWQT(L)) 
          ENDIF
          ! *** MACROALGAE SUBMODEL
          IF(WQI0>0.1 .AND. RMAC(L,K)>0.0)THEN !If there is solar radiation and some macroalgae present in this layer
            IZ=IWQZMAP(L,K)
            WQISM(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPM(IZ)), WQISMIN )  !Optimal light
! *** SOLAR RADIATION AT TOP OF THIS LAYER
            EXPA0=EXP(-WQITOP(L,K)/WQISM(L)) !Macroalgae
! *** UPDATE SOLAR RADIATION AT BOTTOM OF THIS LAYER
            EXPA1=EXP(-WQITOP(L,K)/WQISM(L)*EXP(-WQKESS*DZCHP(L))) !Macroalgae
!*********WQF2IM is the light limitation for macroalgae
		  WQF2IM=EXP(1.0)*WQFD/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0) !Macroalgae light limitation
            XLIMIM(L,K) = XLIMIM(L,K) + WQF2IM !*********XLIMIM is the light limitation for macroalgae
            MACLIM(L,K,2) = WQF2IM !light limitation saved
!*** Maximum macroalgae growth rate modulated by WQF2IM (light limitation}, WQF1NM (nutrient limitation), WQTDGM (temperature limitation), and WQLDF (ecological carrying capacity limitation) 
            WQPM(L)= WQPMM(IMWQZT(L))*WQF2IM*WQF1NM*WQTDGM(IWQT(L))*WQLDF  !Macroalgae growth rate f(I)*h(N)*g(T)*ecological carrying capacity (note that velocity limitation was already considered as a component of nutrient limitation)
            MACLIM(L,K,1) = WQPM(L) !Macroalgae growth rate saved
            !SCJ debug write out!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            if(L==284.and.K==13)then
!                print*,(MACLIM(L,K,NQ),NQ=1,5)!!!!!!!!!!!!!!!!!!!!
!                print*,'Biomass:',WQV(L,K,IDNOTRVA),WQV(L,K,14),WQV(L,19,14)!!!!!!!!!!!
!            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ENDIF  
        ENDDO
C  
C END HORIZONTAL LOOP FOR ALGAE PARMETERS  
C  
        XMRM = 0.0  
        DO L=2,LA  
          IZ=IWQZMAP(L,K)  
          WQOBTOT(L) = WQV(L,K,1)+WQV(L,K,2)+WQV(L,K,3)  
          WQKRPC(L) = (WQKRC + WQKRCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))  
          WQKLPC(L) = (WQKLC + WQKLCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))  
          IF(RMAC(L,K)>0.0) !If macroalgae are present in this layer
     &      XMRM = WQKDCALM(IZ) * WQV(L,K,IDNOTRVA) !Check if it makes sense to multiply by RMAC
C  
C M. MORTON 08/28/99: ADDED SPATIALLY VARIABLE DOC HYDROLYSIS RATE WQKDC  
C    TO ACHIEVE BETTER CONTROL IN SYSTEMS WITH A COMBINATION OF FRESHWAT  
C    STREAMS AND TIDAL RIVERS WITH DIFFERENT CHARACTERISTICS.  
C  
          WQKDOC=(WQKDC(IZ)+WQKDCALG*WQOBTOT(L) + XMRM)*WQTDMNL(IWQT(L))  
          O2WQ(L) = MAX(WQV(L,K,19), 0.0)  
          WQTT1 = WQKDOC / (WQKHORDO + O2WQ(L) + 1.E-18)  
          WQKHR(L) = WQTT1 * O2WQ(L)  
          WQDENIT(L)=WQTT1*WQAANOX*RNO3WQ(L)/(WQKHDNN+RNO3WQ(L) + 1.E-18)  
        ENDDO  
C  
C 7-10 PHOSPHORUS  
C
        ! *** HYDROLYSIS  
        DO L=2,LA  
          WQAPC(L)=1.0/(WQCP1PRM+WQCP2PRM*EXP(-WQCP3PRM*PO4DWQ(L)))  
          WQKHP = (WQKHPC+WQKHPD+WQKHPG) / 3.0  
          WQTT1 = WQKHP / (WQKHP+PO4DWQ(L) + 1.E-18) * WQOBTOT(L)  
          WQKRPP(L) = (WQKRP + WQKRPALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** RPOP--> PO4
          WQKLPP(L) = (WQKLP + WQKLPALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** LPOP--> DOP
          WQKDOP(L) = (WQKDP + WQKDPALG*WQTT1) * WQTDMNL(IWQT(L))  ! *** DOP --> PO4
        ENDDO
        ! *** PHOSPHATE SETTLING   
        DO L=2,LA  
          IF(IWQSRP==1)THEN  
            WQTTM = WQKPO4P*WQTAMP(L,K)  
            WQH10(L) = - WQWSSET(L,1) * WQTTM / (1.0+WQTTM)  
            IF(K/=KC)THEN  
              WQTTM = WQKPO4P*WQTAMP(L,K+1)  
              WQT10(L) = WQWSSET(L,2) * WQTTM / (1.0+WQTTM)  
            ENDIF  
          ELSEIF(IWQSRP==2)THEN  
            WQTTS = WQKPO4P*SEDT(L,K)  
            WQH10(L) = - WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
            IF(K/=KC)THEN  
              WQTTS = WQKPO4P*SEDT(L,K)  
              WQT10(L) = WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
            ENDIF  
          ELSE  
            WQH10(L) = 0.0  
            WQT10(L) = 0.0  
          ENDIF 
          WQH10(L) = WQH10(L)*DTWQO2 
        ENDDO  
C  
C 11-15 NITROGEN  
C  
        ! *** HYDROLYSIS  
        DO L=2,LA  
          WQKHN = (WQKHNC+WQKHND+WQKHNG) / 3.0  
          WQTT1 = WQKHN / (WQKHN+RNH4NO3(L) + 1.E-18) * WQOBTOT(L)  
          WQKRPN(L) = (WQKRN + WQKRNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** RPON-->NH3
          WQKLPN(L) = (WQKLN + WQKLNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** LON -->DON
          WQKDON(L) = (WQKDN + WQKDNALG*WQTT1) * WQTDMNL(IWQT(L))  ! *** DON -->NH3
        ENDDO  
        DO L=2,LA  
          IF(RNH4NO3(L)==0.0)THEN  
            WQPNC(L)=0.0  
            WQPND(L)=0.0  
            WQPNG(L)=0.0  
            WQPNM(L)=0.0  
          ELSE  
            WQTTC = RNH4WQ(L)/(WQKHNC+RNO3WQ(L) + 1.E-18)  
            WQTTD = RNH4WQ(L)/(WQKHND+RNO3WQ(L) + 1.E-18)  
            WQTTG = RNH4WQ(L)/(WQKHNG+RNO3WQ(L) + 1.E-18)  
            WQTTM = RNH4WQ(L)/(WQKHNM+RNO3WQ(L) + 1.E-18)  
            WQPNC(L) = (RNO3WQ(L)/(WQKHNC+RNH4WQ(L) + 1.E-18)  
     &          + WQKHNC/(RNH4NO3(L) + 1.E-18)) * WQTTC  
            WQPND(L) = (RNO3WQ(L)/(WQKHND+RNH4WQ(L) + 1.E-18)  
     &          + WQKHND/(RNH4NO3(L) + 1.E-18)) * WQTTD  
            WQPNG(L) = (RNO3WQ(L)/(WQKHNG+RNH4WQ(L) + 1.E-18)  
     &          + WQKHNG/(RNH4NO3(L) + 1.E-18)) * WQTTG  
            WQPNM(L) = (RNO3WQ(L)/(WQKHNM+RNH4WQ(L) + 1.E-18)  
     &          + WQKHNM/(RNH4NO3(L) + 1.E-18)) * WQTTM  
          ENDIF  
          WQNIT(L) = O2WQ(L) * WQTDNIT(IWQT(L)) /  
     &        ( (WQKHNDO+O2WQ(L)) * (WQKHNN+RNH4WQ(L)) + 1.E-18)  
        ENDDO  
        IF(IWQSI==1)THEN  
          DO L=2,LA  
            IF(IWQSRP==1)THEN  
              WQTTM = WQKSAP*WQTAMP(L,K)  
              WQN17(L) = - WQWSSET(L,1) * WQTTM / (1.0+WQTTM)  
              IF(K/=KC)THEN  
                WQTTM = WQKSAP*WQTAMP(L,K+1)  
                WQT17(L) = WQWSSET(L,2) * WQTTM / (1.0+WQTTM)  
              ENDIF  
            ELSEIF(IWQSRP==2)THEN  
              WQTTS = WQKSAP*SEDT(L,K)  
              WQN17(L) = - WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
              IF(K/=KC)THEN  
                WQTTS = WQKSAP*SEDT(L,K+1)  
                WQT17(L) = WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
              ENDIF  
            ELSE  
              WQN17(L) = 0.0  
              WQT17(L) = 0.0  
            ENDIF  
          ENDDO  
          WQN17(L) = WQN17(L)*DTWQO2 
        ENDIF  
C  
        IF(ISTRWQ(22)==1)THEN
          PPCDO=-3.45 !PARTIAL PRES OF CO2 IN 10^ppcdo ATM; TEMPORARILY DECLARED HERE. SHOULD BE READ IN FROM INPUT FILE
          DO L=2,LA  
            IZ=IWQZMAP(L,K)  
            WQO18(L)= -DTWQO2*WQKCOD(IWQT(L),IZ)*O2WQ(L)/
     &                 (WQKHCOD(IZ) + O2WQ(L) + 1.E-18)
C  
! *** DO Saturation, MOD BY TT, SEE CHAPRA (1997) PG. 3 
!          TVAL1=1./(TWQ(L)+273.15)  
!          TVAL2=TVAL1*TVAL1  
!          TVAL3=TVAL1*TVAL2  
!          TVAL4=TVAL2*TVAL2
!          RLNSAT1=-139.3441+(1.575701E+5*TVAL1)-(6.642308E+7*TVAL2)  
!     &        +(1.2438E+10*TVAL3)-(8.621949E+11*TVAL4)  
!          RLNSAT2=RLNSAT1-SWQ(L)*( 1.7674E-2-(1.0754E+1*TVAL1)  
!     &        +(2.1407E+3*TVAL2) )  
!          WQDOS(L) = EXP(RLNSAT2)  
!          XDOSAT(L,K) = XDOSAT(L,K) + WQDOS(L)*DTWQ*DZCHP(L)   
! *** DO Saturation, Modified by SCJ, see Garcia and Gordon, Limnology and Oceanography 37(6), 1992, Eqn. 8 and Table 1
            TVAL1=LOG((298.15-TWQ(L))/(273.15+TWQ(L)))
            TVAL2=TVAL1*TVAL1
            TVAL3=TVAL1*TVAL2
            TVAL4=TVAL1*TVAL3
            TVAL5=TVAL1*TVAL4
            RLNSAT1=5.80818+3.20684*TVAL1+4.11890*TVAL2+4.93845*TVAL3
     &             +1.01567*TVAL4+1.41575*TVAL5
            RLNSAT2=SWQ(L)*(-7.01211E-3-7.25958E-3*TVAL1-7.93334E-3*TVAL2
     &             -5.54491E-3*TVAL3)-1.32412E-7*SWQ(L)*SWQ(L)
            WQDOS(L) = EXP(RLNSAT1+RLNSAT2)*32E-3 !32E-3 approximately converts micromol/L to mg/L or g/m^3
            XDOSAT(L,K) = XDOSAT(L,K) + WQDOS(L)*DTWQ*DZCHP(L)
          
!************* CO2 parameters
	      CDOSATIDX(L) = -2385.73/(TWQ(L) + 273.15) -	!VB COMPUTING THE pK FOR SAT CONC OF CO2; K - HENRY'S CONST
     &	                       0.0152642 * (TWQ(L) + 273.15) + 14.0184
!          K * MOL WT OF CO2 * PARTIAL PRES OF CO2 IN ATM
	      WQCDOS(L) = 10.**(-CDOSATIDX(L)+PPCDO) * (44.* 1000.) !VB EVALUATING CONC OF CO2 IN G/M^3 
!************* CO2 parameters
          ! *** Compute Reaeration
            IF(K==KC)THEN 
              WINDREA = WINDST(L)  
            ! DO NOT ALLOW WIND SPEEDS ABOVE 11 M/SEC IN THE FOLLOWING EQUATION
              WQWREA=0.728*SQRT(WINDREA)+(0.0372*WINDREA-0.317)*WINDREA  
C  
              LE=LEAST(L)
              LN=LNC(L)
              IF(IWQKA(IZ)==0)THEN  
                WQVREA = WQKRO(IZ)  
                WQWREA = 0.0  
              ELSEIF(IWQKA(IZ)==1)THEN  
                WQVREA = WQKRO(IZ)  
              ELSEIF(IWQKA(IZ)==2)THEN  
                UMRM = 0.5*(U(L,K)+U(LE,K))  
                VMRM = 0.5*(V(L,K)+V(LN,K))  
                XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
              ! *** WQKRO = 3.933 TYPICALLY  
                WQVREA = WQKRO(IZ) * XMRM**0.5 / HP(L)**0.5  
              ELSEIF(IWQKA(IZ)==3)THEN  
                UMRM = MAX(U(L,K), U(LE,K))  
                VMRM = MAX(V(L,K), V(LN,K))  
                XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
              ! *** WQKRO = 5.32 TYPICALLY  
                WQVREA = WQKRO(IZ) * XMRM**0.67 / HP(L)**1.85  
              ELSEIF(IWQKA(IZ)==4)THEN  
            ! *** MODIFIED OWENS AND GIBBS REAERATION EQUATION:  
            ! *** NOTE: NORMALIZED TO A DEPTH OF 1.0 FT, I.E., THIS EQUATION GIVES THE  
            ! ***       SAME REAERATION AS OWENS & GIBBS AT 1.0 FT DEPTH; AT HIGHER  
            ! ***       DEPTHS IT GIVES LARGER REAERATION THAN OWENS & GIBBS.  
                UMRM = MAX(U(L,K), U(LE,K))  
                VMRM = MAX(V(L,K), V(LN,K))  
                XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
                YMRM = HP(L)*3.0*(1.0 - HP(L)/(HP(L)+0.1524))  
            ! *** WQKRO = 5.32 TYPICALLY  
                WQVREA = WQKRO(IZ) * XMRM**0.67 / YMRM**1.85  
              ELSEIF(IWQKA(IZ)== 5)THEN  
                UMRM = MAX(U(L,K), U(LEAST(L),K))  
                VMRM = MAX(V(L,K), V(LNC(L),K))  
                XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
                WQVREA = 3.7*XMRM  
              ENDIF  

            ! *** NOW COMBINE REAERATION DUE TO WATER VELOCITY AND WIND STRESS
              WQVREA = WQVREA * REAC(IZ)  
              WQWREA = WQWREA * REAC(IZ)  
              WQP19(L) = - (WQVREA + WQWREA) * DZWQ(L)* WQTDKR(IWQT(L),IZ)  
              WQKRDOS(L) = -WQP19(L)*WQDOS(L)
              WQP22(L) = WQP19(L)*((32./44.)**0.25) 		!VB Kr FOR CO2 ANALOGOUS TO WQP19 ; 44 = MOL WT OF CO2
              WQKRCDOS(L) = -WQP22(L) * WQCDOS(L)		!VB EVALUATING Kr*SAT CONC OF CO2
            ELSE  
              WQP19(L) = 0.0  
              WQP22(L) = 0.0							!VB Kr FOR CO2 IS ZERO FOR CELLS NOT AT THE SURFACE
            ENDIF  
          ENDDO  
          IF(IWQSRP==1)THEN
            WQR20(2:LA) = WQWPSL(2:LA,K,20)*VOLWQ(2:LA)  
     &          + (WQV(2:LA,K,20) - WQTAMP(2:LA,K)) * WQWSSET(2:LA,1)  
            IF(K==1)THEN
              DO L=2,LA  
                IF(LMASKDRY(L))THEN
                  WQR20(L) = WQR20(L)  
     &            + WQTDTAM(IWQT(L))*DZWQ(L)/(WQKHBMF+O2WQ(L) + 1.E-18)
                ENDIF
              ENDDO
            ENDIF
  
            IF(K/=KC)THEN
              WQR20(2:LA) = WQR20(2:LA)  
     &          + (WQV(2:LA,K+1,20) - WQTAMP(2:LA,K+1)) * WQWSSET(2:LA,2)  
            ELSE     ! K==KC
              WQR20(2:LA)=WQR20(2:LA)+(WQWDSL(2:LA,KC,20)+WQATML(2:LA,KC,20))*VOLWQ(2:LA)
            ENDIF 
          ENDIF
        ENDIF
C  
!      WQA1Cmax=0.0;WQA1Cmin=0.0
        DO M=1,MCOUNT !Macroalgae
          L=IJLMAC(M,3)
          IF(RMAC(L,K)>0.0)THEN
            WQA1C = (WQPM(L) - WQBMM(L) - WQPRM(L)-WQWSM*DZWQ(L))*DTWQO2 !RMAC factor
!           WQA1Cmax=max(WQA1Cmax,WQA1C);WQA1Cmin=min(WQA1Cmin,WQA1C)
            WQVA1C = 1.0 / (1.0 - WQA1C)  
            WQV(L,K,IDNOTRVA)=(WQV(L,K,IDNOTRVA)+WQA1C*WQV(L,K,IDNOTRVA))*WQVA1C !*SMAC(L)  !Macroalgae growth equation
            WQV(L,K,IDNOTRVA) = MAX(WQV(L,K,IDNOTRVA),WQMCMIN) !*SMAC(L)  !Note the lower bound put on macroalgae from Bmin in C44 of WQ3DWC.INP
            WQO(L,IDNOTRVA) = WQVO(L,K,IDNOTRVA)+WQV(L,K,IDNOTRVA)  
          ENDIF
        ENDDO
C
C******************************************************************************
C ***
C *** NOW COMPUTE KINETICS FOR EACH CONSTITUENT
C
C ****  PARAM 01  CHC - cyanobacteria
C
        IF(ISTRWQ(1).EQ.1)THEN  
          DO L=2,LA  
            ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
            WQA1C=(WQPC(L)-WQBMC(L)-WQPRC(L)-WQBCSET(L,1))*DTWQO2 !production per unit time multiplied by half time step
            WQKK(L) = 1.0 / (1.0 - WQA1C)
            ! ***   PT_SRC_LOADS    VOLUME  
            WQR1C = WQWPSL(L,K,1) * VOLWQ(L)  !point source load rate multiplied by inverse cell volume  g/m^3/t
            WQRR(L) = WQV(L,K,1) + DTWQ*WQR1C + WQA1C*WQV(L,K,1)   !transported biomass conc. (CALWQC) + point source load rate * time step + growth rate * previous biomass conc.
          ENDDO  

          IF(K.NE.KC)THEN  
            ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQBCSET(2:LA,2)*WQO(2:LA,1) !biomass conc. + DtX(1/t)* biomass conc.
          ELSE  !Surface layer: K.EQ.KC
            DO L=2,LA  
              ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
              WQR1C = (WQWDSL(L,KC,1)+WQATML(L,KC,1))*VOLWQ(L)  !atmospheric loading mass per time / cell volume
              WQRR(L) = WQRR(L) + DTWQ*WQR1C  !biomass conc. + Dt*loading rate per unit volume
            ENDDO
          ENDIF  
          WQV(2:LA,K,1)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.0-SCB(2:LA))*WQV(2:LA,K,1)  !boundary condition implementation. CAN THIS BE REPLACED BY WQRR(2:LA)*WQKK(2:LA) THROUGHOUT THESE CALCULATIONS?
          WQO(2:LA,1)= WQVO(2:LA,K,1)+WQV(2:LA,K,1) !depth totaled (column sum) biomass conc = old biomass conc in cell + biomass conc from this iteration
        ENDIF  
C
C ****  PARAM 02  CHD - diatom algae
C
        IF(ISTRWQ(2).EQ.1)THEN  
          DO L=2,LA  
            ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
            WQA2D=(WQPD(L)-WQBMD(L)-WQPRD(L)-WQBDSET(L,1))*DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQA2D)  
            ! ***   PT_SRC_LOADS    VOLUME  
            WQR2D = WQWPSL(L,K,2) * VOLWQ(L)  
            WQRR(L) = WQV(L,K,2) + DTWQ*WQR2D + WQA2D*WQV(L,K,2)  
          ENDDO  

          IF(K.NE.KC)THEN  
            ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQBDSET(2:LA,2)*WQO(L,2)
          ELSE  !Surface layer: K.EQ.KC
            DO L=2,LA  
              ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
              WQR2D = (WQWDSL(L,KC,2)+WQATML(L,KC,2))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR2D  
            ENDDO
          ENDIF  
          WQV(2:LA,K,2)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,2)  
          WQO(2:LA,2)=WQVO(2:LA,K,2)+WQV(2:LA,K,2)
        ENDIF  
C
C ****  PARAM 03  CHG - green algae
C
        IF(ISTRWQ(3).EQ.1)THEN
          DO L=2,LA  
            ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
            WQA3G=(WQPG(L)-WQBMG(L)-WQPRG(L)-WQBGSET(L,1))*DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQA3G)  
            ! ***   PT_SRC_LOADS    VOLUME  
            WQR3G = WQWPSL(L,K,3) * VOLWQ(L) 
            ! ***                   External      Internal 
            WQRR(L) = WQV(L,K,3) + DTWQ*WQR3G + WQA3G*WQV(L,K,3)  
          ENDDO  
          IF(K.NE.KC)THEN
            ! *** Add the Algae settled in from the cell above  
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQBGSET(2:LA,2)*WQO(2:LA,3)
          ELSE  !Surface layer: K.EQ.KC
            DO L=2,LA  
              ! ***    ATM DRY DEP   ATM WET DEP     VOLUME  
              WQR3G = (WQWDSL(L,KC,3)+WQATML(L,KC,3))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR3G  
            ENDDO
          ENDIF  
          WQV(2:LA,K,3)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,3)  
          WQO(2:LA,3)=WQVO(2:LA,K,3)+WQV(2:LA,K,3)
        ENDIF  
C
C ****  PARAM 04  ROC - refractory particulate organic carbon
C
        IF(ISTRWQ(4).EQ.1)THEN  
          DO L=2,LA  
            ! ***     HYDROLYSIS    SETTLING  
            WQB4 = -( WQKRPC(L) + WQRPSET(L,1))*DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQB4)  
            ! ***  ALGAE PREDATION SOURCE OF RPOC  
            WQA4 = WQFCRP * (WQPRC(L)*WQO(L,1) + WQPRD(L)*WQO(L,2) + WQPRG(L)*WQO(L,3))  
            IF(RMAC(L,K)>0.0) !If macroalgae is present in this layer
     &        WQA4 = WQA4 + WQFCRPM*WQPRM(L)*WQVO(L,K,IDNOTRVA) !RMAC factor?
            ! ***  PT_SRC_LOADS    VOLUME  
            WQR4 = WQWPSL(L,K,4) * VOLWQ(L)  
            WQRR(L) = WQV(L,K,4) + DTWQ*WQR4 + DTWQO2*WQA4 + WQB4*WQV(L,K,4)  
          ENDDO
  
          IF(K.NE.KC)THEN  
            ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQRPSET(2:LA,2)*WQO(2:LA,4)  
          ELSE  !Surface layer: K.EQ.KC
            DO L=2,LA  
              ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
              WQR4 = (WQWDSL(L,KC,4)+WQATML(L,KC,4))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR4  
            ENDDO
          ENDIF  
          WQV(2:LA,K,4)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,4)  
          WQO(2:LA,4)=WQVO(2:LA,K,4)+WQV(2:LA,K,4)
        ENDIF  
C
C ****  PARAM 05  LOC - labile particulate organic carbon
C
        IF(ISTRWQ(5).EQ.1)THEN  
          DO L=2,LA  
            ! ***     HYDROLYSIS    SETTLING  
            WQC5 = - (WQKLPC(L)  + WQLPSET(L,1))*DTWQO2 
            WQKK(L) = 1.0 / (1.0 - WQC5)    
            WQA5 = WQFCLP * (WQPRC(L)*WQO(L,1) + WQPRD(L)*WQO(L,2) + WQPRG(L)*WQO(L,3)) !Predation
            IF(RMAC(L,K)>0.0) !If macroalgae is present in this layer
     &        WQA5 = WQA5 + WQFCLPM * WQPRM(L)*WQVO(L,K,IDNOTRVA) !RMAC factor?
            ! ***  PT_SRC_LOADS    VOLUME  
            WQR5 = WQWPSL(L,K,5) * VOLWQ(L)
  
            WQRR(L) = WQV(L,K,5) + DTWQ*WQR5 + DTWQO2*WQA5 + WQC5*WQV(L,K,5)    ! *** PMC
          ENDDO  
          IF(K.NE.KC)THEN  
            ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQLPSET(2:LA,2)*WQO(2:LA,5)  
          ELSE  !Surface layer: K.EQ.KC
            DO L=2,LA  
              ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
              WQR5 = (WQWDSL(L,K,5)+WQATML(L,KC,5))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR5  
            ENDDO
          ENDIF  
          WQV(2:LA,K,5)=SCB(L)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,5)  
          WQO(2:LA,5)=WQVO(2:LA,K,5)+WQV(2:LA,K,5)
        ENDIF  
C
C ****  PARAM 06  DOC - dissolved organic carbon
C
        IF(ISTRWQ(6).EQ.1)THEN  
          DO L=2,LA  
            ! ***    RESPIRATION  DENITRIFICATION
            WQD6 = - ( WQKHR(L) +   WQDENIT(L)) *DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQD6)  
            WQA6C=WQFCDC + CFCDCWQ*( WQKHRC/(WQKHRC+O2WQ(L)+ 1.E-18) )  
            WQA6D=WQFCDD + CFCDDWQ*( WQKHRD/(WQKHRD+O2WQ(L)+ 1.E-18) )  
            WQA6G=WQFCDG + CFCDGWQ*( WQKHRG/(WQKHRG+O2WQ(L)+ 1.E-18) )  
            WQA6 = ( WQA6C*WQBMC(L) + WQFCDP*WQPRC(L) )*WQO(L,1)  
     &           + ( WQA6D*WQBMD(L) + WQFCDP*WQPRD(L) )*WQO(L,2)  
     &           + ( WQA6G*WQBMG(L) + WQFCDP*WQPRG(L) )*WQO(L,3)  
            IF(RMAC(L,K)>0.0)THEN !If macroalgae are present in this layer
              IZ=IWQZMAP(L,K)   !NOTE THE INCONSISTENCY HERE WHERE ZONATION IS CONSIDERED. IF ZONATION IS IMPLEMENTED, THE CODE MUST BE UPDATED ACCORDINGLY
              WQA6M=(WQFCDM+(1.-WQFCDM)*WQKHRM(IZ) / (WQKHRM(IZ) + O2WQ(L) + 1.E-18))*WQBMM(L)  
              WQA6 = WQA6 + (WQA6M + WQFCDPM*WQPRM(L))*WQVO(L,K,IDNOTRVA)  
            ENDIF  
            ! ***  PT_SRC_LOADS    VOLUME  
            WQR6 = WQWPSL(L,K,6) * VOLWQ(L)
            WQRR(L) = WQV(L,K,6) + DTWQ*WQR6 + WQD6*WQV(L,K,6) + DTWQO2*(WQA6 +WQKLPC(L)*WQO(L,5))
          ENDDO

          IF(K.EQ.KC)THEN
            DO L=2,LA  
              ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
              WQR6 = (WQWDSL(L,K,6)+WQATML(L,KC,6))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR6  
            ENDDO
          ENDIF
          WQV(2:LA,K,6)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,6)
          WQO(2:LA,6)=WQVO(2:LA,K,6)+WQV(2:LA,K,6)
        ENDIF  
C
C ****  PARAM 07  ROP - refractory particulate organic phosphorus
C
        IF(ISTRWQ(7).EQ.1)THEN  
          DO L=2,LA  
            WQE7 = - (WQKRPP(L)+WQRPSET(L,1))*DTWQO2 
            WQKK(L) = 1.0 / (1.0 - WQE7)  
            WQA7C = (WQFPRC*WQBMC(L) + WQFPRP*WQPRC(L)) * WQO(L,1)  
            WQA7D = (WQFPRD*WQBMD(L) + WQFPRP*WQPRD(L)) * WQO(L,2)  
            WQA7G = (WQFPRG*WQBMG(L) + WQFPRP*WQPRG(L)) * WQO(L,3)  
            WQA7 = (WQA7C+WQA7D+WQA7G) * WQAPC(L)  
            IF(RMAC(L,K)>0.0) !If macroalgae is present in this layer
     &        WQA7 = WQA7 + (WQFPRM*WQBMM(L) + WQFPRPM*WQPRM(L)) * WQVO(L,K,IDNOTRVA)* WQAPC(L)*WQAPCM  
            ! ***  PT_SRC_LOADS    VOLUME  
            WQR7 = WQWPSL(L,K,7) * VOLWQ(L)  
            WQRR(L) = WQV(L,K,7) + DTWQ*WQR7 + DTWQO2*WQA7 + WQE7*WQV(L,K,7)   
          ENDDO  
          IF(K.NE.KC)THEN  
            ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQRPSET(2:LA,2)*WQO(2:LA,7)
          ELSE
            DO L=2,LA  
              ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
              WQR7 = (WQWDSL(L,K,7)+WQATML(L,KC,7))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR7  
            ENDDO
          ENDIF
          WQV(2:LA,K,7)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,7)  
          WQO(2:LA,7)=WQVO(2:LA,K,7)+WQV(2:LA,K,7)
        ENDIF  
C
C ****  PARAM 08  LOP - labile particulate organic phosphorus
C
        IF(ISTRWQ(8).EQ.1)THEN  
          DO L=2,LA  
            ! ***    HYDROLYSIS  SETTLING
            WQF8 = - (WQKLPP(L)+WQLPSET(L,1))*DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQF8)  
            WQA8C = (WQFPLC*WQBMC(L) + WQFPLP*WQPRC(L)) * WQO(L,1)  
            WQA8D = (WQFPLD*WQBMD(L) + WQFPLP*WQPRD(L)) * WQO(L,2)  
            WQA8G = (WQFPLG*WQBMG(L) + WQFPLP*WQPRG(L)) * WQO(L,3)  
            WQA8 = (WQA8C+WQA8D+WQA8G) * WQAPC(L)  
            IF(RMAC(L,K)>0.0) !If macroalgae is present in this layer
     &        WQA8 = WQA8 + (WQFPLM*WQBMM(L) + WQFPLPM*WQPRM(L)) * WQVO(L,K,IDNOTRVA) * WQAPC(L) * WQAPCM
            ! ***  PT_SRC_LOADS    VOLUME  
            WQR8 = WQWPSL(L,K,8) * VOLWQ(L)  
            WQRR(L) = WQV(L,K,8) + DTWQ*WQR8 + DTWQO2*WQA8 + WQF8*WQV(L,K,8)      
          ENDDO  
          IF(K.NE.KC)THEN
            ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQLPSET(2:LA,2)*WQO(2:LA,8)
          ELSE
            DO L=2,LA  
              ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
              WQR8 = (WQWDSL(L,K,8)+WQATML(L,KC,8))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR8  
            ENDDO
          ENDIF  
          WQV(2:LA,K,8)=SCB(2:LA)*(WQRR(L)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,8)  
          WQO(2:LA,8)=WQVO(2:LA,K,8)+WQV(2:LA,K,8)
        ENDIF  
C
C ****  PARAM 09  DOP - dissolved organic phosphorus
C
        IF(ISTRWQ(9).EQ.1)THEN  
          DO L=2,LA
            WQF9 = - DTWQO2*WQKDOP(L)  
            WQKK(L) = 1.0 / (1.0 - WQF9)
            WQA9C = (WQFPDC*WQBMC(L) + WQFPDP*WQPRC(L)) * WQO(L,1)  
            WQA9D = (WQFPDD*WQBMD(L) + WQFPDP*WQPRD(L)) * WQO(L,2)  
            WQA9G = (WQFPDG*WQBMG(L) + WQFPDP*WQPRG(L)) * WQO(L,3)  
            WQA9 = (WQA9C+WQA9D+WQA9G) * WQAPC(L)  
            IF(RMAC(L,K)>0.0) !If macroalgae is present in this layer
     &        WQA9 = WQA9 + (WQFPDM*WQBMM(L) + WQFPDPM*WQPRM(L)) * WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
            ! ***  PT_SRC_LOADS    VOLUME  
            WQR9 = WQWPSL(L,K,9) * VOLWQ(L)  
            WQRR(L) = WQV(L,K,9) + DTWQ*WQR9 + WQF9*WQV(L,K,9) + DTWQO2*(WQA9 + WQKLPP(L)*WQO(L,8) )
          ENDDO

          IF(K.EQ.KC)THEN
            DO L=2,LA  
              ! ***    ATM DRY DEP    ATM WET DEP    VOLUME  
              WQR9 = (WQWDSL(L,KC,9)+WQATML(L,KC,9))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR9  
            ENDDO
          ENDIF
          WQV(2:LA,K,9)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,9)  
        ENDIF  
C
C ****  PARAM 10  P4D - total phosphate
C
        IF(ISTRWQ(10).EQ.1)THEN  
          DO L=2,LA  
            WQA10C=(WQFPIC*WQBMC(L)+WQFPIP*WQPRC(L)-WQPC(L))*WQO(L,1)  
            WQA10D=(WQFPID*WQBMD(L)+WQFPIP*WQPRD(L)-WQPD(L))*WQO(L,2)  
            WQA10G=(WQFPIG*WQBMG(L)+WQFPIP*WQPRG(L)-WQPG(L))*WQO(L,3)
            WQKK(L) = (WQA10C+WQA10D+WQA10G) * WQAPC(L)
            IF(RMAC(L,K)>0.0) !If macroalgae is present in this layer
     &        WQKK(L) = WQKK(L)+(WQFPIM*WQBMM(L)+WQFPIP*WQPRM(L)-WQPM(L))*WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM  !RMAC factor?
          ! ***      PT_SRC_LOADS        VOLUME  
            WQRR(L) = WQWPSL(L,K,10) * VOLWQ(L)  
          ENDDO

          IF(K.EQ.1)THEN  
            DO L=2,LA 
              IF(LMASKDRY(L))THEN !NOTE INCONSISTENCY WITH WETTING/DRYING. IF WETTING/DRYING IS ACTIVE, ALL OF THESE CALCULATIONS SHOULD BE UPDATED ACCORDINGLY
                WQRR(L) = WQRR(L) + WQBFPO4D(L)*DZWQ(L) ! *** Add in Benthic Flux
              ENDIF
            ENDDO  
          ENDIF
          IF(K.EQ.KC)THEN
              ! ***      ATM DRY DEP    ATM WET DEP     VOLUME
!              WQR10 = (WQWDSL(L,KC,10)+WQATML(L,KC,10))*VOLWQ(L)  
!              WQRR(L) = WQRR(L) + WQR10
            DO L=2,LA
              IF(LMASKDRY(L))THEN
              ! ***                  ATM DRY DEP     ATM WET DEP      VOLUME (THESE ARE WQR10)
                WQRR(L) = WQRR(L) + (WQWDSL(L,KC,10)+WQATML(L,KC,10))*VOLWQ(L)
              ENDIF
            ENDDO
          ENDIF
          WQRR(2:LA) = WQV(2:LA,K,10) + DTWQ*WQRR(2:LA) + WQH10(2:LA)*WQV(2:LA,K,10)
     &                       + DTWQO2*(WQKK(2:LA) + WQKRPP(2:LA)*WQO(2:LA,7) + WQKDOP(2:LA)*WQO(2:LA,9))
          ! *** Add in settling from above
          IF(K.NE.KC)WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQT10(2:LA)*WQO(2:LA,10)
          
          DO L=2,LA  
            WQKKL = 1.0 / (1.0 - WQH10(L)) 
            WQV(L,K,10)=SCB(L)*(WQRR(L)*WQKKL)+(1.-SCB(L))*WQV(L,K,10)  
            WQO(L,10)=WQVO(L,K,10)+WQV(L,K,10)
          ENDDO  
        ENDIF  
C
C ****  PARAM 11  RON - refractory particulate organic nitrogen
C
        IF(ISTRWQ(11).EQ.1)THEN  
          DO L=2,LA  
            ! ***     HYDROLYSIS     SETTLING
            WQI11 = - (WQKRPN(L) + WQRPSET(L,1))*DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQI11)  
            WQA11C=(WQFNRC*WQBMC(L)+WQFNRP*WQPRC(L))*WQANCC*WQO(L,1)  
            WQA11D=(WQFNRD*WQBMD(L)+WQFNRP*WQPRD(L))*WQANCD*WQO(L,2)  
            WQA11G=(WQFNRG*WQBMG(L)+WQFNRP*WQPRG(L))*WQANCG*WQO(L,3)  
            WQA11 = WQA11C+WQA11D+WQA11G  
            IF(RMAC(L,K)>0.0) !If macroalgae is in this layer
     &        WQA11 = WQA11 + (WQFNRM*WQBMM(L)+WQFNRPM*WQPRM(L))*WQANCM*WQVO(L,K,IDNOTRVA)  
            ! ***    PT_SRC_LOADS    VOLUME  
            WQR11 = WQWPSL(L,K,11) * VOLWQ(L)
  
            WQRR(L) = WQV(L,K,11) + DTWQ*WQR11 + DTWQO2*WQA11 
     &              + WQI11*WQV(L,K,11) 
          ENDDO 
 
          IF(K.NE.KC)THEN  
            ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQRPSET(2:LA,2)*WQO(2:LA,11)
          ELSE   ! K.EQ.KC
            DO L=2,LA  
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
              WQR11 = (WQWDSL(L,KC,11)+WQATML(L,KC,11))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR11  
            ENDDO
          ENDIF
          WQV(2:LA,K,11)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,11)  
          WQO(2:LA,11)=WQVO(2:LA,K,11)+WQV(2:LA,K,11)
        ENDIF  
C
C ****  PARAM 12  LON - labile particulate organic nitrogen
C
        IF(ISTRWQ(12).EQ.1)THEN  
          DO L=2,LA  
            ! ***     HYDROLYSIS     SETTLING
            WQJ12 = - (WQKLPN(L)+WQLPSET(L,1))*DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQJ12)  
            WQA12C=(WQFNLC*WQBMC(L)+WQFNLP*WQPRC(L))*WQANCC*WQO(L,1)  
            WQA12D=(WQFNLD*WQBMD(L)+WQFNLP*WQPRD(L))*WQANCD*WQO(L,2)  
            WQA12G=(WQFNLG*WQBMG(L)+WQFNLP*WQPRG(L))*WQANCG*WQO(L,3)  
            WQA12 = WQA12C+WQA12D+WQA12G  
            IF(RMAC(L,K)>0.0) !If macroalgae is in this layer
     &        WQA12 = WQA12 + (WQFNLM*WQBMM(L)+WQFNLPM*WQPRM(L))*WQANCM*WQVO(L,K,IDNOTRVA)  
            ! ***    PT_SRC_LOADS    VOLUME  
            WQR12 = WQWPSL(L,K,12) * VOLWQ(L)  

            WQRR(L) = WQV(L,K,12) + DTWQ*WQR12 + DTWQO2*WQA12 + WQJ12*WQV(L,K,12)
          ENDDO
            
          IF(K.NE.KC)THEN  
            ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQLPSET(2:LA,2)*WQO(2:LA,12)
          ELSE   ! K.EQ.KC
            DO L=2,LA  
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
              WQR12 = (WQWDSL(L,KC,12)+WQATML(L,KC,12))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR12  
            ENDDO
          ENDIF  
          WQV(2:LA,K,12)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,12)  
          WQO(2:LA,12)=WQVO(2:LA,K,12)+WQV(2:LA,K,12)
        ENDIF  
C
C ****  PARAM 13  DON - dissolved organic nitrogen
C
        IF(ISTRWQ(13).EQ.1)THEN  
          DO L=2,LA
            WQF13 = - DTWQO2*WQKDON(L)  
            WQKK(L) = 1.0 / (1.0 - WQF13)
            WQA13C=(WQFNDC*WQBMC(L)+WQFNDP*WQPRC(L))*WQANCC*WQO(L,1)  
            WQA13D=(WQFNDD*WQBMD(L)+WQFNDP*WQPRD(L))*WQANCD*WQO(L,2)  
            WQA13G=(WQFNDG*WQBMG(L)+WQFNDP*WQPRG(L))*WQANCG*WQO(L,3)  
            WQA13 = WQA13C+WQA13D+WQA13G  
            IF(RMAC(L,K)>0.0) !If macroalgae is in this layer
     &        WQA13 =WQA13 + (WQFNDM*WQBMM(L)+WQFNDPM*WQPRM(L))*WQANCM*WQVO(L,K,IDNOTRVA) 
            ! ***    PT_SRC_LOADS    VOLUME  
            WQR13 = WQWPSL(L,K,13) * VOLWQ(L)
  
            WQRR(L) = WQV(L,K,13) + DTWQ*WQR13 + WQF13*WQV(L,K,13) + DTWQO2*(WQA13 + WQKLPN(L)*WQO(L,12))
          ENDDO
  
          IF(K.EQ.KC)THEN
            DO L=2,LA  
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
              WQR13 = (WQWDSL(L,KC,13)+WQATML(L,KC,13))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR13  
            ENDDO
          ENDIF  
          WQV(2:LA,K,13)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,13)  
          WQO(2:LA,13)=WQVO(2:LA,K,13)+WQV(2:LA,K,13)
        ENDIF  
C
C ****  PARAM 14  NHX - ammonia nitrogen
C
        IF(ISTRWQ(14).EQ.1)THEN  
          ! ***            PT_SRC_LOADS    VOLUME  
          WQRR(2:LA) = WQWPSL(2:LA,K,14) * VOLWQ(2:LA)  
          IF(K.EQ.1)THEN  
            DO L=2,LA  
              IF(LMASKDRY(L))THEN 
                WQRR(L) = WQRR(L) + WQBFNH4(L)*DZWQ(L)   ! *** Add in Benthic Flux
              ENDIF
            ENDDO  
          ELSEIF(K.EQ.KC)THEN
            DO L=2,LA  
              IF(LMASKDRY(L))THEN
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
                WQR14 = (WQWDSL(L,KC,14)+WQATML(L,KC,14))*VOLWQ(L)  
                WQRR(L) = WQRR(L) + WQR14  
              ENDIF
            ENDDO
          ENDIF  

          DO L=2,LA
            WQF14 = - DTWQO2*WQNIT(L)  
            WQKK(L) = 1.0 / (1.0 - WQF14) 
            WQA14C=WQFNIC*WQBMC(L)+WQFNIP*WQPRC(L)-WQPNC(L)*WQPC(L)  
            WQA14D=WQFNID*WQBMD(L)+WQFNIP*WQPRD(L)-WQPND(L)*WQPD(L)  
            WQA14G=WQFNIG*WQBMG(L)+WQFNIP*WQPRG(L)-WQPNG(L)*WQPG(L)  
            WQA14 = WQA14C*WQANCC*WQO(L,1) + WQA14D*WQANCD*WQO(L,2) + WQA14G*WQANCG*WQO(L,3)  
            IF(RMAC(L,K)>0.0) !If macroalgae is in this layer
     &         WQA14 = WQA14 + (WQFNIM*WQBMM(L)+WQFNIPM*WQPRM(L) - WQPNM(L)*WQPM(L))*WQANCM*WQVO(L,K,IDNOTRVA)
            WQRR(L) = WQV(L,K,14) + DTWQ*WQRR(L) + WQF14*WQV(L,K,14) + DTWQO2*(WQA14 + WQKRPN(L)*WQO(L,11) + WQKDON(L)*WQO(L,13))        
            WQV(L,K,14)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,14)  
            WQO(L,14)=WQVO(L,K,14)+WQV(L,K,14)
          ENDDO  
        ENDIF  
C
C ****  PARAM 15  NOX - nitrate nitrogen
C
        IF(ISTRWQ(15).EQ.1)THEN  
            ! ***      PT_SRC_LOADS    VOLUME  
          WQRR(2:LA) = WQWPSL(2:LA,K,15) * VOLWQ(2:LA)  
          IF(K.EQ.1)THEN  
            DO L=2,LA  
              IF(LMASKDRY(L))THEN 
                WQRR(L) = WQRR(L) + WQBFNO3(L)*DZWQ(L)   ! *** Add in Benthic Flux
              ENDIF
            ENDDO  
          ELSEIF(K.EQ.KC)THEN
            DO L=2,LA  
              IF(LMASKDRY(L))THEN
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
                WQR15 = (WQWDSL(L,KC,15)+WQATML(L,KC,15))*VOLWQ(L)  
                WQRR(L) = WQRR(L) + WQR15  
              ENDIF
            ENDDO
          ENDIF  
          DO L=2,LA  
            WQA15C = (WQPNC(L)-1.0)*WQPC(L) * WQANCC * WQO(L,1)  
            WQA15D = (WQPND(L)-1.0)*WQPD(L) * WQANCD * WQO(L,2)  
            WQA15G = (WQPNG(L)-1.0)*WQPG(L) * WQANCG * WQO(L,3)  
            WQA15 = WQA15C+WQA15D+WQA15G  
            IF(RMAC(L,K)>0.0) !If macroalgae is in this layer
     &        WQA15 =WQA15 + (WQPNM(L)-1.0)*WQPM(L)*WQANCM*WQVO(L,K,IDNOTRVA)
            WQB15 = WQV(L,K,15) + DTWQ*WQRR(L) + DTWQO2*(WQA15 - WQANDC*WQDENIT(L)*WQO(L,6) + WQNIT(L)*WQO(L,14))

            WQV(L,K,15)=SCB(L)*WQB15 + (1.-SCB(L))*WQV(L,K,15)  
            WQO(L,15)=WQVO(L,K,15)+WQV(L,K,15)
          ENDDO  
        ENDIF  
C
C ****  PARAM 16  SUU - particulate biogenic silica
C
        IF(ISTRWQ(16).EQ.1.AND.IWQSI.EQ.1)THEN  
          DO L=2,LA  
            WQM16 = - (WQKSUA(IWQT(L)) + WQBDSET(L,1)) * DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQM16)  
            WQA16D = (WQFSPD*WQBMD(L) + WQFSPP*WQPRD(L)) * WQASCD * WQO(L,2)  
            ! ***    PT_SRC_LOADS    VOLUME  
            WQR16 = WQWPSL(L,K,16) * VOLWQ(L)

            WQRR(L) = WQV(L,K,16) + DTWQ*WQR16 + DTWQO2*WQA16D + WQM16*WQV(L,K,16)  
          ENDDO  
          IF(K.NE.KC)THEN  
          ! *** Add in settling from above
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQBDSET(2:LA,2)*WQO(2:LA,16)     ! *** PMC
          ELSE
            DO L=2,LA  
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
              WQR16 = (WQWDSL(L,KC,16)+WQATML(L,KC,16))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR16  
            ENDDO
          ENDIF  
          WQV(2:LA,K,16)=SCB(2:LA)*WQRR(2:LA)*WQKK(2:LA) + (1.-SCB(2:LA))*WQV(2:LA,K,16)  
          WQO(2:LA,16)=WQVO(2:LA,K,16)+WQV(2:LA,K,16)
        ENDIF
C
C ****  PARAM 17  SAA - dissolved available silica
C
        IF(ISTRWQ(17).EQ.1.AND.IWQSI.EQ.1)THEN
          WQKK(2:LA) = (WQFSID*WQBMD(2:LA) + WQFSIP*WQPRD(2:LA) - WQPD(2:LA)) * WQASCD * WQO(2:LA,2)  
          ! ***      PT_SRC_LOADS    VOLUME  
          WQRR(2:LA) = WQWPSL(2:LA,K,17) * VOLWQ(2:LA)  
          
          IF(K.EQ.1)THEN  
            DO L=2,LA
              IF(LMASKDRY(L))THEN 
                WQRR(L) = WQRR(L) + WQBFSAD(L)*DZWQ(L)   ! *** Add in Benthic Flux
              ENDIF
            ENDDO  
          ENDIF  
          
          WQRR(2:LA)= WQV(2:LA,K,17) + 
     &      DTWQ*WQRR(2:LA) +WQN17(2:LA)*WQV(2:LA,K,17) + DTWQO2*(WQKK(2:LA) + WQKSUA(IWQT(2:LA))*WQO(2:LA,16))  
          
          IF(K.NE.KC)THEN  
            WQRR(2:LA) = WQRR(2:LA) + DTWQO2*WQT17(2:LA)*WQO(2:LA,17)  
          ELSE
            DO L=2,LA  
              ! ***      ATM DRY DEP     ATM WET DEP    VOLUME  
              WQR17 = (WQWDSL(L,KC,17)+WQATML(L,KC,17))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR17  
            ENDDO
          ENDIF  
          WQKK(2:LA) = 1.0 / (1.0 - WQN17(2:LA))
          WQV(2:LA,K,17)=SCB(2:LA)*WQRR(2:LA)*WQKK(2:LA) + (1.0-SCB(2:LA))*WQV(2:LA,K,17)  
          WQO(2:LA,17)=WQVO(2:LA,K,17)+WQV(2:LA,K,17)
        ENDIF  
C
C ****  PARAM 18  COD - chemical oxygen demand
C
        IF(ISTRWQ(18).EQ.1)THEN  
          WQKK(2:LA) = 1.0 / (1.0 - WQO18(2:LA))  
              ! ***       PT_SRC_LOADS       VOLUME  
          WQRR(2:LA) = WQWPSL(2:LA,K,18) * VOLWQ(2:LA)
          IF(K.EQ.1)THEN  
            DO L=2,LA  
              IF(LMASKDRY(L))THEN 
                WQRR(L) = WQRR(L) + WQBFCOD(L)*DZWQ(L)   ! *** Add in Benthic Flux
              ENDIF
            ENDDO    
          ELSEIF(K.EQ.KC)THEN
            DO L=2,LA  
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
              WQR18 = (WQWDSL(L,KC,18)+WQATML(L,KC,18))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + WQR18  
            ENDDO
          ENDIF  
          WQRR(2:LA)=WQV(2:LA,K,18)+DTWQ*WQRR(2:LA)+WQO18(2:LA)*WQV(2:LA,K,18)  
          WQV(2:LA,K,18)=SCB(2:LA)*(WQRR(2:LA)*WQKK(2:LA))+(1.-SCB(2:LA))*WQV(2:LA,K,18)  
          WQO(2:LA,18) = WQVO(2:LA,K,18)+WQV(2:LA,K,18)  
        ENDIF
C  
C ****  PARAM 19  DOX - dissolved oxygen
C
        ! ***  1) CHC - cyanobacteria
        ! ***  2) CHD - diatom algae
        ! ***  3) CHG - green algae
        ! ***  4) ROC - refractory particulate organic carbon
        ! ***  5) LOC - labile particulate organic carbon
        ! ***  6) DOC - dissolved organic carbon
        ! ***  7) ROP - refractory particulate organic phosphorus
        ! ***  8) LOP - labile particulate organic phosphorus
        ! ***  9) DOP - dissolved organic phosphorus
        ! *** 10) P4D - total phosphate
        ! *** 11) RON - refractory particulate organic nitrogen 
        ! *** 12) LON - labile particulate organic nitrogen
        ! *** 13) DON - dissolved organic nitrogen
        ! *** 14) NHX - ammonia nitrogen
        ! *** 15) NOX - nitrate nitrogen
        ! *** 16) SUU - particulate biogenic silica
        ! *** 17) SAA - dissolved available silica
        ! *** 18) COD - chemical oxygen demand
        ! *** 19) DOX - dissolved oxygen
        ! *** 20) TAM - total active metal
        ! *** 21) FCB - fecal coliform bacteria
        ! *** 22) CO2 - dissolved carbon dioxide
        ! *** 23) MAC - macroalgae

        !C 04/29/99 MRM:  
        !C THE FOLLOWING ARRAYS WERE ADDED TO KEEP TRACK OF THE VARIOUS COMPONENT  
        !C OF DISSOLVED OXYGEN.  THE INSTANTANEOUS VALUES FOR EACH COMPONENT ARE  
        !C SUMMED IN THE ARRAYS AND THEN DUMPED TO THE WQDOCOMP.BIN FILE AT THE  
        !C SAME TIME INTERVAL AS FOR THE WQWCAVG.BIN FILES (I.E., IWQTSDT INTERVA  
        !C USUALLY DAILY AVERAGES).  THE ARRAY DESCRIPTIONS ARE:  
        !C  XDOPSL(L,L) = D.O. COMPONENT FOR POINT SOURCE LOADS
        !C  XDOSOD(L,K) = D.O. COMPONENT FOR SEDIMENT OXYGEN DEMAND  
        !C  XDOKAR(L,K) = D.O. COMPONENT FOR REAERATION  
        !C  XDODOC(L,K) = D.O. COMPONENT FOR DISS. ORG. CARBON DECAY  
        !C  XDONIT(L,K) = D.O. COMPONENT FOR AMMONIA NITRIFICATION  
        !C  XDOCOD(L,K) = D.O. COMPONENT FOR CHEM. OXY. DEMAND OXIDATION  
        !C  XDOPPB(L,K) = D.O. COMPONENT FOR PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
        !C  XDORRB(L,K) = D.O. COMPONENT FOR RESPIRATION OF TOTAL CHLOROPHYLL  
        !C  XDOPPM(L,K) = D.O. COMPONENT FOR PHOTOSYNTHESIS OF MACROALGAE  
        !C  XDORRM(L,K) = D.O. COMPONENT FOR RESPIRATION OF MACROALGAE  
        !C  XDOALL(L,K) = SUM OF THE ABOVE 10 D.O. COMPONENTS  
        !C  NLIM = COUNTER FOR NUMBER OF ITEMS SUMMED IN EACH ARRAY SLOT  
      
        IF(ISTRWQ(19).EQ.1)THEN  
          DO L=2,LA  
            WQRR(L) = WQWPSL(L,K,19) * VOLWQ(L)  
            TMP19=WQRR(L)*DTWQ*DZCHP(L)
            XDOPSL(L,K) = XDOPSL(L,K) + TMP19
            XDOALL(L,K) = XDOALL(L,K) + TMP19  
          ENDDO  

          ! *** Handle Surface Processes
          IF(K.EQ.KC)THEN  
            DO L=2,LA  
              WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L)) 
              ! ***             ATM DRY DEP    ATM WET DEP    VOLUME  
              WQRR(L)=WQRR(L)+(WQWDSL(L,KC,19)+WQATML(L,KC,19))*VOLWQ(L) !VJ - volume should multiply both load terms
              ! *** Reaeration
              !WQRR(L) = WQRR(L) - WQP19(L)*(WQDOS(L) - WQV(L,K,19)) !This is the equation from WQ theory manual
	        ! O2-reaeration changes
              WQRR(L) = WQRR(L) + WQKRDOS(L)  !WQKRDOS=-WQP19*WQDOS SCJ this seems wrong because it allows too much DO if it is produced by a lot of algae
              WQV(L,KC,19)=MIN(WQV(L,KC,19),WQDOS(L)) !This does not allow O2 to exceed saturation at the water surface for high density algae systems
            ENDDO
          ELSE 
            WQKK(2:LA) = 1.0
          ENDIF 

          ! *** Bottom Processes 
          IF(K.EQ.1)THEN  
            DO L=2,LA  
              IF(LMASKDRY(L))THEN 
                TEMFAC=1.065**(TEM(L,1)-20.)  
                WQRR(L) = WQRR(L) + TEMFAC*WQBFO2(L)*DZWQ(L)   ! *** Add in Benthic Flux
                TMP19=TEMFAC*WQBFO2(L)*DTWQ
                XDOSOD(L,K) = XDOSOD(L,K) + TMP19  
                XDOALL(L,K) = XDOALL(L,K) + TMP19
              ENDIF
            ENDDO  
          ENDIF  
C 
          DO L=2,LA  
            DTWQxH = DTWQ*DZCHP(L)  
            DTWQxH2= DTWQO2*DZCHP(L)

            IF(WQI0 .LE. 0.001)THEN  
              WQTTC = 0.0
              WQTTD = 0.0
              WQTTG = 0.0
            ELSE
              WQTTC = (1.3 - 0.3*WQPNC(L)) * WQPC(L) 
              WQTTD = (1.3 - 0.3*WQPND(L)) * WQPD(L)  
              WQTTG = (1.3 - 0.3*WQPNG(L)) * WQPG(L)

              ! *** PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
              TMP19 = WQAOCR*DTWQxH2*
     &          (WQTTC*WQO(L,1)+WQTTD*WQO(L,2)+WQTTG*WQO(L,3))
              XDOPPB(L,K) = XDOPPB(L,K) + TMP19
              XDOALL(L,K) = XDOALL(L,K) + TMP19
            ENDIF

            ! *** RESPIRATION OF TOTAL CHLOROPHYLL - CYANOBACTERIA 
            XMRM = CFCDCWQ*O2WQ(L)*WQBMC(L)/(WQKHRC+O2WQ(L)+ 1.E-18)  
            WQA19C = WQTTC - XMRM
            TMP19  = XMRM*WQO(L,1)*WQAOCR * DTWQxH2  
            XDORRB(L,K) = XDORRB(L,K) - TMP19 
            XDOALL(L,K) = XDOALL(L,K) - TMP19

            ! *** RESPIRATION OF TOTAL CHLOROPHYLL - DIATOMS 
            XMRM = CFCDDWQ*O2WQ(L)*WQBMD(L)/(WQKHRD+O2WQ(L)+ 1.E-18)  
            WQA19D = WQTTD - XMRM
            TMP19  = XMRM*WQO(L,2)*WQAOCR * DTWQxH2
            XDORRB(L,K) = XDORRB(L,K) - TMP19  
            XDOALL(L,K) = XDOALL(L,K) - TMP19

            ! *** RESPIRATION OF TOTAL CHLOROPHYLL -  GREENS
            XMRM = CFCDGWQ*O2WQ(L)*WQBMG(L)/(WQKHRG+O2WQ(L)+ 1.E-18)  
            WQA19G = WQTTG - XMRM
            TMP19  = XMRM*WQO(L,3)*WQAOCR * DTWQxH2   
            XDORRB(L,K) = XDORRB(L,K) - TMP19
            XDOALL(L,K) = XDOALL(L,K) - TMP19

            ! *** TOTAL NET RESPIRATION/PHOTOSYNTHESIS
            WQA19=(WQA19C*WQO(L,1) + WQA19D*WQO(L,2) + WQA19G*WQO(L,3))
     &             * WQAOCR
            ! *** MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED  
            ! ***   TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE  
            IF(RMAC(L,K)>0.0)THEN !If macroalgae is in this layer
              ! *** TRAPEZOIDAL AVERAGE CONCENTRATIONS
              WQO(L,IDNOTRVA)=WQVO(L,K,IDNOTRVA)+WQV(L,K,IDNOTRVA)
              IZ = IWQZMAP(L,K)  
              WQTTM = (1.3 - 0.3*WQPNM(L)) * WQPM(L)  
              XMRM=(1.0-WQFCDM)*O2WQ(L)*WQBMM(L)/(WQKHRM(IZ)+O2WQ(L) + 1.E-18) 
              WQA19A = WQTTM*WQO(L,IDNOTRVA)*WQAOCRPM - XMRM*WQO(L,IDNOTRVA)*WQAOCRRM 
              WQA19 = WQA19 + WQA19A
              TMP19 = WQTTM*WQO(L,IDNOTRVA)*WQAOCRPM * DTWQxH2     
              XDOPPM(L,K) = XDOPPM(L,K) + TMP19
              XDOALL(L,K) = XDOALL(L,K) + TMP19 

              TMP19 = XMRM*WQO(L,IDNOTRVA)*WQAOCRRM * DTWQxH2
              XDORRM(L,K) = XDORRM(L,K) - TMP19
              XDOALL(L,K) = XDOALL(L,K) - TMP19
            ENDIF  

            ! *** O2 Mass Balance
            ! WQA19                         ! *** Total Net Respiration/Photosynthesis
            WQSUM=DTWQ*WQRR(L)              ! *** Sum of Loadings/Demands
            WQRea=WQP19(L)*WQV(L,K,19)      ! *** Reaeration
            WQPOC=WQAOCR*WQKRPC(L)*WQO(L,4) ! *** POC
            WQDOC=WQAOCR*WQKHR(L) *WQO(L,6) ! *** DOC
            WQNH3=WQAONT*WQNIT(L) *WQO(L,14)! *** Ammonia
            WQCOD=WQO18(L)*WQO(L,18)        ! *** COD
            WQRR(L) = WQV(L,K,19) + WQSUM  + WQCOD  + 
     &                DTWQO2*(WQA19 - WQPOC - WQDOC - WQNH3 + WQRea)
c            WQRR(L) = WQVO(L,K,19) + DTWQ*WQRR(L) +   
c     &                DTWQO2*( WQA19 - WQAOCR*WQKRPC(L)*WQVO(L,K,4)  
c     &                               - WQAOCR*WQKHR(L) *WQVO(L,K,6)
c     &                               - WQAONT*WQNIT(L) *WQVO(L,K,14)  
c     &                               + WQP19(L)*WQVO(L,K,19) )
c     &                               + WQO18(L)*WQVO(L,K,18)  
            WQV(L,K,19)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,19)  ! *** remove the SCB computations!
C            WQV(L,K,19)=SCB(L)*WQRR(L)+(1.-SCB(L))*WQV(L,K,19)  ! *** remove the SCB computations!
            ! *** WQV(L,K,19) After this WQV(L,K,19) can not be < 0.
            WQV(L,K,19) = MAX(WQV(L,K,19), 0.0)  
            !WQVO(L,K,19) = WQVO(L,K,19)+WQV(L,K,19)  
C  
            ! *** COMPUTE AND SAVE D.O. DEFICIT:  
            IF(ISMTSB.LT.ISMTSE)THEN  
              WQO(L,19)=WQVO(L,K,19)+WQV(L,K,19)
              XMRM = WQDOS(L) - WQV(L,K,19)  
              XDODEF(L,K) = XDODEF(L,K) + XMRM*DTWQ*DZCHP(L)  
              IF(K.EQ.KC)THEN
                TMP19=WQKRDOS(L)*DTWQ*DZCHP(L) + 
     &                WQP19(L)*WQO(L,19)*DTWQxH2  
                XDOKAR(L,K) = XDOKAR(L,K) + TMP19  
                XDOALL(L,K) = XDOALL(L,K) + TMP19  
              ENDIF

              TMP19=WQAOCR*WQKHR(L)*WQO(L,6)*DTWQxH2  
              XDODOC(L,K)=XDODOC(L,K) - TMP19  
              XDOALL(L,K)=XDOALL(L,K) - TMP19  

              TMP19=WQAONT*WQNIT(L)*WQO(L,14)*DTWQxH2  
              XDONIT(L,K)=XDONIT(L,K) - TMP19
              XDOALL(L,K)=XDOALL(L,K) - TMP19

              TMP19=WQO18(L)*WQO(L,18)*DZCHP(L)
              XDOCOD(L,K)=XDOCOD(L,K) - TMP19  
              XDOALL(L,K)=XDOALL(L,K) - TMP19  

              XDODZ(L,K) = XDODZ(L,K) + DZCHP(L)  
            ENDIF
		! O2-reaeration changes
		! If DO conc. is greater than saturated DO conc., excess DO
		! gets released out, which is calculated as XDOOUT.
		  IF(K==KC)THEN
		    IF(WQV(L,KC,19)>WQDOS(L))THEN
			  XDOOUT(L,KC)=(WQV(L,KC,19)-WQDOS(L))*DTWQ*DZCHP(L)
			  WQV(L,KC,19)=WQDOS(L)
			  XDOALL(L,K)=XDOALL(L,K)-XDOOUT(L,KC)
              ENDIF
		  ENDIF

          ENDDO  
        ENDIF  
C
C ****  PARAM 20  TAM - total active metal
C
        IF(ISTRWQ(20).EQ.1)THEN  
          IF(IWQSRP.EQ.1)THEN  
            DO L=2,LA  
              WQT20 = - DTWQ*WQWSSET(L,1)    ! *** DTWQO2
              WQKK(L) = 1.0 / (1.0 - WQT20)  
              WQRR(L)=WQV(L,K,20)+DTWQ*WQR20(L)+WQT20*WQV(L,K,20)  
            ENDDO  
            IF(K.NE.KC)THEN  
              ! *** Add in settling from above
              DO L=2,LA  
                WQRR(L) = WQRR(L) + DTWQO2*WQWSSET(L,2)*WQO(L,20)
              ENDDO  
            ENDIF  
              WQV(2:LA,K,20)=SCB(2:LA)*( WQRR(2:LA)*WQKK(2:LA) )  
     &            +(1.-SCB(2:LA))*WQV(2:LA,K,20)  
              WQO(2:LA,20)=WQVO(2:LA,K,20)+WQV(2:LA,K,20)
          ENDIF  
        ENDIF  
C  
C ****  PARAM 21  FCB - fecal coliform bacteria
C
        IF(ISTRWQ(21).EQ.1)THEN  
          IF(IWQFCB.EQ.1)THEN  
            DO L=2,LA  
              WQKK(L) = WQTD2FCB(IWQT(L))  
C  
              ! ***      ATM DRY DEP     LOADS        VOLUME  
		    WQR21= (WQWDSL(L,K,21)+WQWPSL(L,K,21))*VOLWQ(L)
C  
              IF(K==KC .AND. LMASKDRY(L))THEN
              ! ***               ATM WET DEP     VOLUME
                WQR21 = WQR21 + (WQATML(L,KC,21)*VOLWQ(L))  
              ENDIF
		    WQRR(L) = WQV(L,K,21)*WQTD1FCB(IWQT(L)) + DTWQ*WQR21
              WQV(L,K,21)=SCB(L)*( WQRR(L)*WQKK(L) )  
     &            +(1.-SCB(L))*WQV(L,K,21)  
              !WQVO(L,K,21) = WQVO(L,K,21)+WQV(L,K,21)  
            ENDDO  
          ENDIF  
        ENDIF
C  
C ****  PARAM 22 DISSOLVED CARBON DIOXIDE
C
!C THE FOLLOWING ARRAYS WERE ADDED TO KEEP TRACK OF THE VARIOUS COMPONENT  
!C OF DISSOLVED CARBON DIOXIDE.  
!C THE ARRAY DESCRIPTIONS ARE:  
!C  XCDOKAR(L,K) = CDO. COMPONENT FOR REAERATION  
!C  XCDODOC(L,K) = CDO. COMPONENT FOR DISS. ORG. CARBON DECAY  
!C  XCDOPPB(L,K) = CDO. COMPONENT FOR PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
!C  XCDORRB(L,K) = CDO. COMPONENT FOR RESPIRATION OF TOTAL CHLOROPHYLL  
!C  XCDOPPM(L,K) = CDO. COMPONENT FOR PHOTOSYNTHESIS OF MACROALGAE  
!C  XCDORRM(L,K) = CDO. COMPONENT FOR RESPIRATION OF MACROALGAE  
!C  XCDOALL(L,K) = SUM OF THE ABOVE 6 CDO. COMPONENTS  

	  IF(ISTRWQ(22).EQ.1)THEN  
          WQRR(2:LA) = WQWPSL(2:LA,K,22) * VOLWQ(2:LA)  
          ! *** Handle Surface Processes
          IF(K.EQ.KC)THEN  
            DO L=2,LA  
              WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP22(L)) 
          ! ***             ATM DRY DEP    ATM WET DEP    VOLUME  
              WQRR(L)=WQRR(L)+(WQWDSL(L,KC,22)+WQATML(L,KC,22))*VOLWQ(L)
          ! *** Reaeration
              WQRR(L) = WQRR(L) + WQKRCDOS(L)  
            ENDDO
          ELSE  
            WQKK(2:LA) = 1.0
          ENDIF 	
          DO L=2,LA  
            DTWQxH = DTWQ*DZCHP(L)  
            DTWQxH2= DTWQO2*DZCHP(L)
            IF(WQI0 .LE. 0.001)THEN  
              WQTTC = 0.0
              WQTTD = 0.0
              WQTTG = 0.0
            ELSE
              WQTTC = (1.3 - 0.3*WQPNC(L)) * WQPC(L) 
              WQTTD = (1.3 - 0.3*WQPND(L)) * WQPD(L)  
              WQTTG = (1.3 - 0.3*WQPNG(L)) * WQPG(L)
        ! *** PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
            ENDIF
        ! *** RESPIRATION OF TOTAL CHLOROPHYLL - CYANOBACTERIA 
            XMRM = CFCDCWQ*O2WQ(L)*WQBMC(L)/(WQKHRC+O2WQ(L)+ 1.E-18)  
            WQA22C = WQTTC - XMRM
        ! *** RESPIRATION OF TOTAL CHLOROPHYLL - DIATOMS 
            XMRM = CFCDDWQ*O2WQ(L)*WQBMD(L)/(WQKHRD+O2WQ(L)+ 1.E-18)  
            WQA22D = WQTTD - XMRM
        ! *** RESPIRATION OF TOTAL CHLOROPHYLL -  GREENS
            XMRM = CFCDGWQ*O2WQ(L)*WQBMG(L)/(WQKHRG+O2WQ(L)+ 1.E-18)  
            WQA22G = WQTTG - XMRM
        ! *** TOTAL NET RESPIRATION/PHOTOSYNTHESIS
            WQA22=3.67*(WQA22C*WQO(L,1)+WQA22D*WQO(L,2)+WQA22G*WQO(L,3))  !VB 3.67 CONVERTS g CARBON TO g CO2
        ! *** MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED  
        ! *** TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE  
        ! *** TRAPEZOIDAL AVERAGE CONCENTRATIONS
        ! *** CO2 Mass Balance
        ! WQA22												! *** Total Net Respiration/Photosynthesis
            WQCDSUM=DTWQ*WQRR(L)							! *** Sum of Loadings/Demands
            WQCDRea=WQP22(L)*WQV(L,K,22)					! *** Reaeration
            WQCDDOC=(WQKHR(L)+WQDENIT(L))*WQO(L,6)*3.67		! *** DOC FROM HYDROLYSIS AND DENITRIFICATION 3.67 CONVERTS G CARBON TO G CO2    

            WQRR(L) = WQV(L,K,22) + WQCDSUM  +  
     &                DTWQO2*(-WQA22 + WQCDDOC + WQCDRea)
            WQV(L,K,22)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,22)  ! *** remove the SCB computations!
           
        ! *** WQV(L,K,22) After this WQV(L,K,22) can not be < 0.
            WQV(L,K,22) = MAX(WQV(L,K,22), 0.0)  
          ENDDO  
        ENDIF  
      ENDDO  ! *** END OF THE KC LOOP
C ----------------------------------------------------------------  
C  
C INCREMENT COUNTER FOR LIMITATION AND XDOXXX DO COMPONENT ARRAYS:  
C  
      IF(ISDYNSTP.EQ.0)THEN  
        TIMTMP=DT*FLOAT(N)+TCON*TBEGIN  
        TIMTMP=TIMTMP/TCTMSR  
      ELSE  
        TIMTMP=TIMESEC/TCTMSR  
      ENDIF  
      TIMESUM3 = TIMESUM3 + TIMTMP  
      NLIM = NLIM + 1  

C PMC - Moved CHLa Computations to the beginning of the WQ Calculations
      IF(IWQSRP.EQ.1)THEN  
        ! *** Sorption Option: TAM
        DO K=1,KC  
          DO L=2,LA  
            O2WQ(L) = MAX(WQV(L,K,19), 0.0)  
            WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ(L)), WQV(L,K,20) )  
            WQTAMP(L,K) = WQV(L,K,20) - WQTAMD  
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*WQTAMP(L,K))  
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*WQTAMP(L,K))  
          ENDDO  
        ENDDO  
      ELSE IF(IWQSRP.EQ.2)THEN
        ! *** Sorption Option: Sediments 
        WQPO4D(2:LA,1:KC) = WQV(2:LA,1:KC,10) / (1.0 + WQKPO4P*SEDT(2:LA,1:KC))
        WQSAD(2:LA,1:KC)  = WQV(2:LA,1:KC,17) / (1.0 + WQKSAP*SEDT(2:LA,1:KC))
      ELSE  
        WQPO4D(2:LA,1:KC) = WQV(2:LA,1:KC,10)
        WQSAD(2:LA,1:KC)  = WQV(2:LA,1:KC,17)
      ENDIF  
C  
C COUPLING TO SEDIMENT MODEL  
C EVALUATE DEP. FLUX USING NEW VALUES CAUSE IMPLICIT SCHEME IS USED IN SPM
C  
      IF(IWQBEN.EQ.1)THEN  
        DO L=2,LA  
          IMWQZ = IWQZMAP(L,1)  
          WQDFBC(L) = SCB(L)*WQWSC(IMWQZ)*WQV(L,1,1)  
          WQDFBD(L) = SCB(L)*WQWSD(IMWQZ)*WQV(L,1,2)  
          WQDFBG(L) = SCB(L)*WQWSG(IMWQZ)*WQV(L,1,3)  
     &        +WQWSM*DZWQ(L)*WQV(L,1,IDNOTRVA)  
          WQDFRC(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,1,4)  
          WQDFLC(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,1,5)  
          WQDFRP(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,1,7)  
          WQDFLP(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,1,8)  
          WQDFRN(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,1,11)  
          WQDFLN(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,1,12)  
          IF(IWQSI.EQ.1) WQDFSI(L) = SCB(L)*WQWSD(IMWQZ)*WQV(L,1,16)  
        ENDDO  
        IF(IWQSRP.EQ.1)THEN  
          DO L=2,LA  
            IMWQZ = IWQZMAP(L,1)  
            WQDFLP(L) = SCB(L)*( WQDFLP(L)  
     &          + WQWSS(IMWQZ)*( WQV(L,1,10)-WQPO4D(L,1) ) )  
            IF(IWQSI.EQ.1) WQDFSI(L) = SCB(L)*( WQDFSI(L)  
     &          + WQWSS(IMWQZ)*( WQV(L,1,17)-WQSAD(L,1) ) )  
          ENDDO  
        ELSE IF(IWQSRP.EQ.2)THEN  
          DO L=2,LA  
            WQDFLP(L) = SCB(L)*( WQDFLP(L)+WSEDO(NS)*( WQV(L,1,10)  
     &          -WQPO4D(L,1) ) )  
            IF(IWQSI.EQ.1) WQDFSI(L) = SCB(L)*( WQDFSI(L)  
     &          + WSEDO(NS)*( WQV(L,1,17)-WQSAD(L,1) ) )  
          ENDDO  
        ENDIF  
      ENDIF  
C  
C DIURNAL DO ANALYSIS  
C  
      IF(NDDOAVG.GE.1.AND.DEBUG)THEN  
        OPEN(1,FILE='DIURNDO.OUT',POSITION='APPEND')  
        NDDOCNT=NDDOCNT+1  
        NSTPTMP=NDDOAVG*NTSPTC/2  
        RMULTMP=1./FLOAT(NSTPTMP)  
        DDOMAX(2:LA,1:KC)=MAX(DDOMAX(2:LA,1:KC),WQV(2:LA,1:KC,19))
        DDOMIN(2:LA,1:KC)=MIN(DDOMIN(2:LA,1:KC),WQV(2:LA,1:KC,19))
        IF(NDDOCNT.EQ.NSTPTMP)THEN  
          NDDOCNT=0  
          IF(ISDYNSTP.EQ.0)THEN  
            TIME=DT*FLOAT(N)+TCON*TBEGIN  
            TIME=TIME/TCON  
          ELSE  
            TIME=TIMESEC/TCON  
          ENDIF  
          WRITE(1,1111)N,TIME  
          DO L=2,LA  
            WRITE(1,1112)IL(L),JL(L),(DDOMIN(L,K),K=1,KC),  
     &          (DDOMAX(L,K),K=1,KC)  
          ENDDO  
          DDOMAX(2:LA,1:KC)=-1.0E6
          DDOMIN(2:LA,1:KC)= 1.0E6
        ENDIF  
        CLOSE(1)  
      ENDIF
      ! *** LIGHT EXTINCTION ANALYSIS  
      IF(NDLTAVG.GE.1)THEN  
        OPEN(1,FILE='LIGHT.OUT',POSITION='APPEND')  
        NDLTCNT=NDLTCNT+1  
        NSTPTMP=NDLTAVG*NTSPTC/2  
        RMULTMP=1./FLOAT(NSTPTMP)  
        DO K=1,KC  
          DO L=2,LA  
            RLIGHT1=WQKEB(IMWQZT(L))+WQKETSS*SEDT(L,K)  
            XMRM = WQKECHL*WQCHL(L,K)  
            IF(WQKECHL .LT. 0.0)THEN  
              XMRM = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)  
            ENDIF  
            RLIGHT2 = XMRM  
            RLIGHTT(L,K)=RLIGHTT(L,K)+RLIGHT1  
            RLIGHTC(L,K)=RLIGHTC(L,K)+RLIGHT1+RLIGHT2  
          ENDDO  
        ENDDO  
        IF(NDLTCNT.EQ.NSTPTMP)THEN  
          NDLTCNT=0  
          IF(ISDYNSTP.EQ.0)THEN  
            TIME=DT*FLOAT(N)+TCON*TBEGIN  
            TIME=TIME/TCON  
          ELSE  
            TIME=TIMESEC/TCON  
          ENDIF  
          RLIGHTT(2:LA,1:KC)=RMULTMP*RLIGHTT(2:LA,1:KC)
          RLIGHTC(2:LA,1:KC)=RMULTMP*RLIGHTC(2:LA,1:KC)
          WRITE(1,1111)N,TIME  
          DO L=2,LA  
            WRITE(1,1113)IL(L),JL(L),(RLIGHTT(L,K),K=1,KC),  
     &          (RLIGHTC(L,K),K=1,KC)  
          ENDDO  
          RLIGHTT(2:LA,1:KC)=0.0
          RLIGHTC(2:LA,1:KC)=0.0
        ENDIF  
        CLOSE(1)  
      ENDIF
      !PRINT*,WQV(LIJ(17,2),13,15),WQV(LIJ(17,3),13,15),WQV(LIJ(17,10),13,15),WQV(LIJ(17,11),13,15)

 1111 FORMAT(I12,F10.4)  
 1112 FORMAT(2I5,12F7.2)  
 1113 FORMAT(2I5,12E12.4)  
 1414 FORMAT(I12,11E12.4)  
      RETURN  
      END
