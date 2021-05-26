
      SUBROUTINE HDMT2T  
C  
C **  SUBROUTINE HDMT2T EXECUTES THE FULL HYDRODYNAMIC AND MASS  
C **   TRANSPORT TIME INTERGATION USING A TWO TIME LEVEL SCHEME  
C  
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION  
C  
C----------------------------------------------------------------------C  
C  
C CHANGE RECORD  
C DATE MODIFIED     BY                 DATE APPROVED    BY  
C  
C 05/01/2002        John Hamrick       05/01/2002       John Hamrick  
C  modified calls to calbal and budget subroutines  
C  added calls to bal2t1, bal2t4, bal2t5  
C 05/02/2002        John Hamrick       05/01/2002       John Hamrick  
C  modified calculation of cell center bed stress (stored as QQ(l,0))  
C  for cells have source/sinks  
C 09-22-2004        Paul M. Craig  
C  Merged DS and TT versions with the 06-04-2004 TT code  
C----------------------------------------------------------------------C  
C  
C**********************************************************************C  
C  
 
! This is the preferred option for checking for nan using ISNAN
! For older versions of gfortran (< 5.1),  replace with isnan
!      use ieee_arithmetic

      USE GLOBAL  
      USE DRIFTER
      USE PARALLEL_MPI
      USE WINDWAVE ,ONLY:WINDWAVEINIT,WINDWAVETUR
      USE OMP_LIB
#ifdef key_ncdf
      USE IOM
#endif

#ifdef key_mpi
      USE mpi
#endif

#ifdef key_omp
      REAL(8) ostart,oend,tomp,lstart,lend,lomp
#else
      REAL(4) ostart,oend
#endif
      REAL  TMP,TCHKROFF, VFW2
      
      INTEGER::I1,I2
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::ISSBCP  
      LOGICAL BTEST, LTEST, CELL_INSIDE_DOMAIN, KINVARA_MC

      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WCOREW  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::WCORNS  

      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::LCORNER  
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::LCORNWE  
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::LCORNSN  

! arrays for computing tidal range
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::MXR  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::MNR  

      REAL, DIMENSION(2):: DYNAMIC_TIME_STEP_ARRAY_LOCAL, DYNAMIC_TIME_STEP_ARRAY_REDUCE
      REAL SAL_CELL
C
!      INTERFACE TO FUNCTION KBHIT  
!     &    [C,ALIAS:'__kbhit']  
!     &    ()  
!      LOGICAL KBHIT*1  
!      END  
!      INTERFACE TO FUNCTION GETCH  
!     &    [C,ALIAS:'__getch']  
!     &    ()  
!      INTEGER GETCH*1  
!      END  
C
      KINVARA_MC = .FALSE. ! This activates a check on freshwater mass balance and prints to file
      KINSALE_DILUTION = .FALSE. ! This activates computation of dilution coefficient
      COMP_RESIDENCE = .FALSE.
      IF (KINVARA_MC) THEN
         OPEN(777, FILE='CheckFreshwaterMassBal.txt', status='unknown')
         VFW2 = 0.
         IFLPRINT = 0
         VSW_BOUNDARY = 0
      END IF
      IF(.NOT.ALLOCATED(WCOREW))THEN
        ALLOCATE(WCOREW(LCM))
        ALLOCATE(WCORNS(LCM))
        ALLOCATE(LCORNER(LCM))
        ALLOCATE(LCORNWE(LCM))
        ALLOCATE(LCORNSN(LCM))
        ALLOCATE(MXR(LCM))
        ALLOCATE(MNR(LCM))

        ! *** ALLOCATE LOCAL ARRAYS
        WCOREW=0.0 
        WCORNS=0.0
        LCORNER=0
        LCORNWE=0
        LCORNSN=0
        MXR = 0.
        MNR = 0.
      ENDIF
C  
      TTMP=SECNDS(SECNDS_ZERO)  
      ICALLTP=0
C  
      ISTL=2  
      FOURDPI=4./PI  
      ISTL=2   
      IS2TL=1
  
C**********************************************************************C  
C  
C **  SET FLAGS FOR CORNER CELL BED STRESS CORRECTIONS  
C  
C *** DSLLC BEGIN BLOCK
      IF(ISCORTBC.GE.1) THEN  
C  
C **  SET FLAG FOR CELLS HAVING VOLUME SOURCE OR SINKS  
C  
        IF(.NOT.ALLOCATED(ISSBCP))ALLOCATE(ISSBCP(LCM))  
        DO L=1,LC  
          ISSBCP(L)=0  
        ENDDO  
C  
        DO L=2,LA  
          IF(RSSBCE(L).GT.1.5)ISSBCP(L)=1  
          IF(RSSBCW(L).GT.1.5)ISSBCP(L)=1  
          IF(RSSBCN(L).GT.1.5)ISSBCP(L)=1  
          IF(RSSBCS(L).GT.1.5)ISSBCP(L)=1  
        ENDDO  
C
      ENDIF
C  
      DO L=2,LA  
        WCOREST(L)=1.  
        WCORWST(L)=1.  
        WCORNTH(L)=1.  
        WCORSTH(L)=1.  
      ENDDO  
      ! *** DSLLC
C  
C**********************************************************************C  
C  
C **  REINITIALIZE VARIABLES  
C  
      DO L=2,LA  
        H1P(L)=HP(L)  
        H1U(L)=HU(L)  
        H1UI(L)=HUI(L)  
        H1V(L)=HV(L)  
        H1VI(L)=HVI(L)  
        UHDY1E(L)=UHDYE(L)  
        VHDX1E(L)=VHDXE(L)  
      ENDDO  
C  
      DO K=1,KC  
        DO L=2,LA  
          U1(L,K)=U(L,K)  
          V1(L,K)=V(L,K)  
          UHDY1(L,K)=UHDY(L,K)  
          VHDX1(L,K)=VHDX(L,K)  
        ENDDO  
      ENDDO  
C  
C**********************************************************************C  
C  
C **  INITIALIZE COURANT NUMBER DIAGNOSTICS  
C  
      DO K=1,KC  
        DO L=2,LA  
          CFLUUU(L,K)=0.  
          CFLVVV(L,K)=0.  
          CFLWWW(L,K)=0.  
          CFLCAC(L,K)=0.  
        ENDDO  
      ENDDO  
C  
C**********************************************************************C  
C  
      ILOGC=0  
C  
C**********************************************************************C  
C  
C **  CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING  
C **  CALCULATE VELOCITY GRADIENTS  
C  
C----------------------------------------------------------------------C  
C  
      DO L=2,LA  
        LN=LNC(L)  
        LS=LSC(L)
        LE=LEAST(L)
        LW=LWEST(L)
        LNW=LNWC(L)  
        LSE=LSEC(L)  
        LSW=LSWC(L)  
  
        UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))  
     &             +HP(L )*(U(LE ,1)+U(L ,1)))*HVI(L)  
        U1V(L)=0.25*(H1P(LS)*(U1(LSE,1)+U1(LS,1))  
     &              +H1P(L )*(U1(LE ,1)+U1(L ,1)))*H1VI(L)  
        VU(L)=0.25*(HP(LW)*(V(LNW,1)+V(LW,1))  
     &             +HP(L )*(V(LN ,1)+V(L ,1)))*HUI(L)  
        V1U(L)=0.25*(H1P(LW)*(V1(LNW,1)+V1(LW,1))  
     &              +H1P(L )*(V1(LN ,1)+V1(L ,1)))*H1UI(L)  
        ! *** DSLLC END BLOCK
      ENDDO  

C  
C**********************************************************************C  
C  
C **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS  
C  
      IF(ISWAVE.EQ.1) CALL WAVEBL    
      IF(ISWAVE.EQ.2) CALL WAVESXY  
      IF(ISWAVE.EQ.3.AND.NWSER > 0) THEN
        CALL WINDWAVEINIT  
        CALL WINDWAVETUR   !DHC FIRST CALL
      ENDIF
C  
C**********************************************************************C  
C  
C **  FIRST CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS  
C  
      DTDYN=DT  ! *** PMC - FOR INITIALIZATION
      CALL CALTBXY(ISTL,IS2TL)  
C  
C**********************************************************************C  
C  
C **  CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES  
C  
      IF(ISHDMF.GE.1) CALL CALHDMF  
C  
C**********************************************************************C  
C  
C **  CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N  
C  
C----------------------------------------------------------------------C  
C  
      N=-1  
      ! This call to CALTSXY used primarily to initialize wind
      ! sheltering coefficents. Not needed for DT coupling
      CALL CALTSXY
C  
C**********************************************************************C  
C  
C **  SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS  
C  
      CALL CALTBXY(ISTL,IS2TL)  
C  
C**********************************************************************C  
C  
C **  SET BOTTOM AND SURFACE STRESSES  
C  
C----------------------------------------------------------------------C  
C  
      DO L=2,LA  
        TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)  
     &      +U(L,1)*U(L,1)))*U(L,1)  
        TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)  
     &      +V(L,1)*V(L,1)))*V(L,1)  
      ENDDO  
C  
      N=0  
      CALL CALTSXY  
C  
C----------------------------------------------------------------------C  
C  
C **  SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES  
C  
      DO L=2,LA  
        HDFUFX(L)=1.  
        HDFUFY(L)=1.  
        HDFUF(L)=1.  
      ENDDO  
C  
      IF(ISBSDFUF.GE.1)THEN  
        HDFUFM=1.E-12  
C  
        DO L=2,LA  
          LS=LSC(L)  
          HDFUFX(L)=HDFUFM+G*SUB(L)*HU(L)*(BELV(LWEST(L))-BELV(L))*DXIU(L)  
          HDFUFY(L)=HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )-BELV(L))*DYIV(L)  
        ENDDO  
C  
        DO L=2,LA  
          HDFUFX(L)=TBX(L)/HDFUFX(L)  
          HDFUFY(L)=TBY(L)/HDFUFY(L)  
        ENDDO  
C  
        DO L=2,LA  
          HDFUFX(L)=MAX(HDFUFX(L),-1.0)  
          HDFUFY(L)=MAX(HDFUFY(L),-1.0)  
        ENDDO  
C  
        DO L=2,LA  
          HDFUFX(L)=MIN(HDFUFX(L),1.0)  
          HDFUFY(L)=MIN(HDFUFY(L),1.0)  
        ENDDO  
C  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED  
C  
C----------------------------------------------------------------------C  
C  
      IF(ISWAVE.EQ.0)THEN  
C  
C----------------------------------------------------------------------c  
C  
        IF(ISCORTBC.EQ.0) THEN  
C  
          DO L=2,LA  
            TVAR3S(L)=TSY(LNC(L))  
            TVAR3W(L)=TSX(LEAST(L))  
            TVAR3E(L)=TBX(LEAST(L)   )  
            TVAR3N(L)=TBY(LNC(L))  
          ENDDO  
C  
          DO L=2,LA  
            QQ(L,0 )=0.5*CTURB2*SQRT(  
     &          (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2  
     &          +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2)  
            TAUBSED(L)=QQ(L,0 )/CTURB2  
            QQ(L,KC)=0.5*CTURB2*SQRT(  
     &          (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2  
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)  
            QQSQR(L,0)=SQRT(QQ(L,0))  ! *** DSLLC
          ENDDO  
C  
        ENDIF  
C  
C----------------------------------------------------------------------c  
C  
        IF(ISCORTBC.GE.1) THEN  
C  
          IF(DEBUG)THEN
            IF(ISCORTBCD.GE.1)THEN  
              OPEN(1,FILE='ADJSTRESSE.OUT')  
              CLOSE(1,STATUS='DELETE')  
            ENDIF  
C  
            OPEN(1,FILE='TBCORINIT.OUT')  
          ENDIF
C  
          DO L=2,LA  
            IF(ISSBCP(L).EQ.0)THEN  
              IF(SUB(LEAST(L)).LT.0.5) WCOREST(L)=FSCORTBCV(L)  
              IF(SUB(L).LT.0.5) WCORWST(L)=FSCORTBCV(L)  
              IF(SVB(LNC(L)).LT.0.5) WCORNTH(L)=FSCORTBCV(L)  
              IF(SVB(L).LT.0.5) WCORSTH(L)=FSCORTBCV(L)  
            ENDIF  
          ENDDO  
C  
          DO L=2,LA  
            WCOREW(L)=1./(WCOREST(L)+WCORWST(L))  
            WCORNS(L)=1./(WCORNTH(L)+WCORSTH(L))  
          ENDDO  
C  
          DO L=2,LA  
            WCOREST(L)=WCOREST(L)*WCOREW(L)  
            WCORWST(L)=WCORWST(L)*WCOREW(L)  
            WCORNTH(L)=WCORNTH(L)*WCORNS(L)  
            WCORSTH(L)=WCORSTH(L)*WCORNS(L)  
          ENDDO  
C  
          DO L=2,LA  
            TVAR3S(L)=TSY(LNC(L))  
            TVAR3W(L)=TSX(LEAST(L))  
            TVAR3E(L)=TBX(LEAST(L)   )  
            TVAR3N(L)=TBY(LNC(L))  
          ENDDO  
C  
          DO L=2,LA  
            QQ(L,0 )=CTURB2*SQRT(  
     &          (RSSBCE(L)*WCOREST(L)*TVAR3E(L)  
     &          +RSSBCW(L)*WCORWST(L)*TBX(L))**2  
     &          +(RSSBCN(L)*WCORNTH(L)*TVAR3N(L)  
     &          +RSSBCS(L)*WCORSTH(L)*TBY(L))**2)  
            QQ(L,KC)=0.5*CTURB2*SQRT(  
     &          (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2  
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)  
            QQSQR(L,0)=SQRT(QQ(L,0))  ! *** DSLLC
          ENDDO  
C  
          DO L=2,LA  
            TAUBSED(L)=QQ(L,0 )/CTURB2  
            TAUBSND(L)=QQ(L,0 )/CTURB2  
          ENDDO  
  
          IF(DEBUG)THEN
            DO L=2,LA  
              IF(WCORSTH(L).LT.0.49.OR.WCORSTH(L).GT.0.51)THEN  
                IF(WCORWST(L).LT.0.49.OR.WCORWST(L).GT.0.51)THEN  
                  WRITE(1,3678)IL(L),JL(L),WCORWST(L),WCOREST(L),  
     &              WCORSTH(L),WCORNTH(L)  
                ENDIF  
              ENDIF  
            ENDDO  
C  
            CLOSE(1)  
          ENDIF
C  
        ENDIF  
C  
C----------------------------------------------------------------------c  
C  
      ENDIF  
C  
C     ENDIF  
C  
C**********************************************************************C  
C  
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED  
C  
C----------------------------------------------------------------------C  
C  
      IF(ISWAVE.GE.1)THEN  
C  
        DO L=2,LA  
          TVAR3S(L)=TSY(LNC(L))  
          TVAR3W(L)=TSX(LEAST(L))  
          TVAR3E(L)=TBX(LEAST(L)   )  
          TVAR3N(L)=TBY(LNC(L))  
        ENDDO  
C  
        DO L=2,LA  
          TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2  
     &            +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2  
          TAUBC=0.5*SQRT(TAUBC2)  
          UTMP=0.5*STCUV(L)*(U(LEAST(L),1)+U(L,1))+1.E-12  
          VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))  
          CURANG=ATAN2(VTMP,UTMP)  
          TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))  
     &        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))  
          TAUB2=MAX(TAUB2,0.)  
          QQ(L,0 )=CTURB2*SQRT(TAUB2)  
          QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2  
     &        +(TVAR3S(L)+TSY(L))**2)  
          QQSQR(L,0)=SQRT(QQ(L,0))  ! *** DSLLC
        ENDDO  
C  
      ENDIF  
C  
C     ENDIF  
C  
C**********************************************************************C  
C  
C **   SET SWITCHES FOR TWO TIME LEVEL INTEGRATION  
C  
      ISTL=2  
      IS2TL=1  
      DELT=DT  
      DELTD2=DT/2.  
      DZDDELT=DZ/DELT  
C  
C**********************************************************************C  
C  
C **  BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT  
C **  CALCULATION  
C  
C **  SET CYCLE COUNTER AND CALL TIMER  
C  
      NTIMER=0  
      ISSREST=0
      NRESTO=ISRESTO*NTSPTC  
      N=0  
C  
C *** EE BEGIN BLOCK  
C **  INITIALZE & RECORD TIME  
C  
      TIMEDAY=TCON*TBEGIN/86400.  
      NYRS = FLOOR(TIMEDAY/365.)  
      TCHKROFF=(NYRS + 1) *365.
C     CALL TIMELOG(0,TIMEDAY)  
      IF(ISDYNSTP.GT.0)THEN 
        ! *** ALLOW FOR SEDIMENT RAMPUP 
        ISAVESEDDT=ISEDDT   
        ISEDDT=1
      ENDIF
C  
C *** EE END BLOCK  
C  
      NTIMER=1  
      NINCRMT=1  
      NLOOP=0
C  
C----------------------------------------------------------------------C  
C  
        if (partid == MASTER_TASK) CALL CPU_TIME(ostart) 
 1001 CONTINUE  
      IF(N.GE.NTS)GO TO 1000  
            


C  

      IF(ISDYNSTP.EQ.0)THEN  
        N=N+1  
        ETIMESEC=DT*FLOAT(N)  
        ETIMEDAY=DT*FLOAT(N)/86400.  
        TIMESEC=(DT*FLOAT(N)+TCON*TBEGIN)  
        TIMEDAY=(DT*FLOAT(N)+TCON*TBEGIN)/86400.  
      ELSE  
        NLOOP=NLOOP+1  
        IF(NLOOP.GT.ITRMADJ)THEN
          NLOOP=0
          ISEDDT=ISAVESEDDT ! *** PMC-ALLOW FOR SEDIMENT RAMPUP ALSO
          IF(IDRYTBP.EQ.0)THEN  
            CALL CALSTEP  
          ELSE  
            CALL CALSTEPD  
          ENDIF  
          !! UPDATE THE DYNAMCIALLY COMPUTED TIMESTEP FROM EACH PROCESSOR
          !! CALL A REDUCE ACROSS PROCESSORS, IF MPI PARALLEL
#ifdef key_mpi
          ! We need to syncronize timestep and nincrmt across processors
          DYNAMIC_TIME_STEP_ARRAY_LOCAL(1) = NINCRMT
          DYNAMIC_TIME_STEP_ARRAY_LOCAL(2) = DTDYN
          CALL MPI_ALLREDUCE(DYNAMIC_TIME_STEP_ARRAY_LOCAL,DYNAMIC_TIME_STEP_ARRAY_REDUCE,
     &                      2,MPI_REAL8,MPI_MIN,EFDC_COMM,IERR1)
          NINCRMT = DYNAMIC_TIME_STEP_ARRAY_REDUCE(1)
          DTDYN = DYNAMIC_TIME_STEP_ARRAY_REDUCE(2)  ! Set dynamic timestep equal to minimum
#endif
        ELSE
          DTDYN=DT
          NINCRMT=1
        ENDIF
        DELT=DTDYN  
        DELTD2=DTDYN/2.  
        DZDDELT=DZ/DTDYN  
        N=N+NINCRMT
        ETIMESEC=DT*FLOAT(N)  
        ETIMEDAY=(DT*FLOAT(N))/86400.  
        TIMESEC=(DT*FLOAT(N)+TCON*TBEGIN)  
        TIMEDAY=(DT*FLOAT(N)+TCON*TBEGIN)/86400.  
      ENDIF  
C  
C PMC      IF(ILOGC.EQ.NTSMMT)THEN  
      IF(ILOGC.EQ.NTSPTC)THEN  
        CLOSE(8,STATUS='DELETE') 
        OPEN(8,FILE='EFDCLOG'//ans(partid2)//'.OUT',STATUS='UNKNOWN')   
        IF(DEBUG)THEN 
          IF(ISDRY.GT.0)THEN  
            OPEN(1,FILE='DRYWET.LOG',STATUS='UNKNOWN')  
            CLOSE(1,STATUS='DELETE')  
          ENDIF  
          IF(ISCFL.EQ.1)THEN
            OPEN(1,FILE='CFL.OUT',STATUS='UNKNOWN')  
            CLOSE(1,STATUS='DELETE')  
          ENDIF  
        ENDIF
        ILOGC=0  
      ENDIF  
C  
      IF(ISDYNSTP.EQ.0)THEN  
        ILOGC=ILOGC+1  
      ELSE  
        ILOGC=ILOGC+NINCRMT  
      ENDIF  
C  
C *** DSLLC BEGIN BLOCK
      IF(N.LE.NLTS)THEN
        SNLT=0.
      ELSEIF(N.GT.NLTS.AND.N.LE.NTTS)THEN  
        NTMP1=N-NLTS  
        NTMP2=NTTS-NLTS+1  
        SNLT=FLOAT(NTMP1)/FLOAT(NTMP2)  
      ELSE
        SNLT=1.  
      ENDIF  
      
      ! *** TURN OFF WIND SHELTERING FOR ICE CONDITIONS (TO BE REPLACED AFTER FULL ICE SUBMODEL ADDED)
       IF(WINTER_END > WINTER_START)THEN
        IF(TIMEDAY > WINTER_START)THEN
          IF(WINDSTKA_SAVE(1)==0.)THEN
            ! *** TOGGLE OFF THE WIND SHELTERING COEFFICIENTS
            DO L=2,LA
              WINDSTKA_SAVE(L)=WINDSTKA(L)
              WINDSTKA(L)=0.
            ENDDO 
            WINDSTKA_SAVE(1) = 1.
! Fearghal: introduce first generation ice model parameterisation here
! parametrization based on Wang (2010) Development of the Great Lakes Ice Circulation Model (GLIM) - Journal Great Lakes Research
          REVC_SAVE=REVC; RCHC_SAVE=RCHC
          REVC=0.0; RCHC=0.0
          DO L=2,LA  
            CLEVAP(L)=0.001*ABS(REVC)
            CCNHTT(L)=0.001*ABS(RCHC)  
          ENDDO  
          SWRATNF_SAVE=SWRATNF; SWRATNS_SAVE=SWRATNS
          SWRATNF=0.28*SWRATNF; SWRATNS=0.28*SWRATNS
          ENDIF
          IF(TIMEDAY > WINTER_END)THEN
            ! *** TOGGLE ON THE WIND SHELTERING COEFFICIENTS
            DO L=2,LA
              WINDSTKA(L) = WINDSTKA_SAVE(L)
            ENDDO 
            REVC=REVC_SAVE; RCHC=RCHC_SAVE
            DO L=2,LA  
              CLEVAP(L)=0.001*ABS(REVC)
              CCNHTT(L)=0.001*ABS(RCHC)  
            ENDDO  
            SWRATNF=SWRATNF_SAVE; SWRATNS=SWRATNS_SAVE
            WINDSTKA_SAVE(1) = 0.
            WINTER_START = WINTER_START+365.
            WINTER_END = WINTER_END+365.
          ENDIF
        ENDIF  
      ENDIF
C *** DSLLC END BLOCK
C  
      IF(N.LE.NTSVB)THEN  
        GP=GPO*(FLOAT(N)/FLOAT(NTSVB))  
      ELSE  
        GP=GPO  
      ENDIF 
C  
C----------------------------------------------------------------------C  
C  
C **  INITIALIZE TWO-TIME LEVEL BALANCES  
C  
      IF(IS2TIM.GE.1) THEN  
        IF(ISBAL.GE.1)THEN  
          CALL BAL2T1  
        ENDIF  
      ENDIF  
C  
C----------------------------------------------------------------------C  
C  
C **  REENTER HERE FOR TWO TIME LEVEL LOOP
C  
  500 CONTINUE  
C  
C**********************************************************************C  
C  
C **  CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)  
C  
      T1TMP=SECNDS(SECNDS_ZERO)  
      IF(KC.GT.1)THEN  
        IF(ISQQ.EQ.1)THEN 
          IF(ISTOPT(0).EQ.0)CALL CALAVBOLD (ISTL)  
          IF(ISTOPT(0).GE.1)CALL CALAVB (ISTL)
        ELSEIF(ISQQ.EQ.2)THEN
          CALL CALAVB2 (ISTL)  
        ENDIF
      ENDIF  
      TAVB=TAVB+SECNDS(T1TMP)  
C  
C**********************************************************************C  
C  
C **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS  
C  
      IF(ISWAVE.EQ.1) CALL WAVEBL  
      IF(ISWAVE.EQ.2) CALL WAVESXY  
      IF(ISWAVE.EQ.3.AND.NWSER > 0) CALL WINDWAVETUR   !DHC NEXT CALL

C  
C**********************************************************************C  
C  
C **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND UPDATE NEW WIND   
C **  STRESSES  *** DSLLC MOVED 
C  
C----------------------------------------------------------------------C  
C
      CALL CALTSXY  

 
C**********************************************************************C  
C  
C **  CALCULATE EXPLICIT MOMENTUM EQUATION TERMS  
C  
      T1TMP=SECNDS(SECNDS_ZERO)  
      IF(IS2TIM.EQ.1) CALL CALEXP2T  
      IF(IS2TIM.EQ.2) CALL CALIMP2T  
      TCEXP=TCEXP+SECNDS(T1TMP)
C  
C**********************************************************************C  
C  
C **  UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS, CONCENTRATIONS,  
C **  VEGETATION CHARACTERISTICS AND SURFACE ELEVATIONS  
C 
      CALL CALCSER (ISTL)  
      CALL CALVEGSER (ISTL)  
      CALL CALQVS (ISTL)  
      PSERT(0)=0.  
      IF(NPSER.GE.1) CALL CALPSER (ISTL)  
C
C**********************************************************************C
C
C **  WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
C     SETUP CALL
C
C----------------------------------------------------------------------C
C
      IF(ISWSEDA.GT.0.OR.ISUVDA.GT.0) CALL PUVDASM(ISTL,1)
C
C  
C**********************************************************************C  
C  
C **  SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE  
C  
      T1TMP=SECNDS(SECNDS_ZERO)  
      IF(ISCHAN.EQ.0.AND.ISDRY.EQ.0) CALL CALPUV2T 
      IF(ISCHAN.GE.1.OR.ISDRY.GE.1) CALL CALPUV2C  
      TPUV=TPUV+SECNDS(T1TMP)  
C  
C**********************************************************************C  
C  
C **  WRITE DIAGNOSTICS  
C  
C----------------------------------------------------------------------C  
C  
C **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE  
C **  SUPPORTED ON OTHER SYSTEMS.  
C  
      IF(ISLOG.GE.1)THEN  
        WRITE(8,17)N,ITER,RSQ,CFMAX,AVMAX,ABMIN,ABMAX,ABMIN  
      ENDIF  
C  
   17 FORMAT('  N,ITER,RSQ,CFMAX,AVMAX,AVMIN,ABMAX,ABMIN',  
     &    I7,I5,2E12.4,4(1X,F8.4))  
C  
      ERRMAX=MAX(ERRMAX,ERR)  
      ERRMIN=MIN(ERRMIN,ERR)  
      ITRMAX=MAX(ITRMAX,ITER)  
      IRRMIN=MIN(ITRMIN,ITER)  
C  
C**********************************************************************C  
C  
C **  ADVANCE INTERNAL VARIABLES  
C  
C----------------------------------------------------------------------C  
C  
      DO K=1,KC  
        DO L=2,LA  
          UHDY2(L,K)=UHDY1(L,K)  
          UHDY1(L,K)=UHDY(L,K)  
          VHDX2(L,K)=VHDX1(L,K)  
          VHDX1(L,K)=VHDX(L,K)  
          U2(L,K)=U1(L,K)  
          V2(L,K)=V1(L,K)  
          U1(L,K)=U(L,K)  
          V1(L,K)=V(L,K)  
          W2(L,K)=W1(L,K)  
          W1(L,K)=W(L,K)  
        ENDDO  
      ENDDO  
C  
C**********************************************************************C  
C  
C **  SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W  
C  
C----------------------------------------------------------------------C  
C  

      T1TMP=SECNDS(SECNDS_ZERO)  
      IF(KC.GT.1)THEN
        CALL CALUVW (ISTL,IS2TL) 
      ELSE  
        DO L=2,LA  
          UHDY(L,1)=UHDYE(L)  
          U(L,1)=UHDYE(L)*HUI(L)*DYIU(L)  
          VHDX(L,1)=VHDXE(L)  
          V(L,1)=VHDXE(L)*HVI(L)*DXIV(L)  
          W(L,1)=0.  
        ENDDO  
        CALL CALUVW (ISTL,IS2TL)  
      ENDIF  
      TUVW=TUVW+SECNDS(T1TMP)  
C**********************************************************************C
C
C **  WATER SURFACE ELEVATION AND VELOCITY DATA ASSIMILATION
C     DIAGNOSTIC CALL
C
C----------------------------------------------------------------------C
C
      IF(ISWSEDA.GT.0.OR.ISUVDA.GT.0) CALL PUVDASM(ISTL,2)
C
C  
C**********************************************************************C  
C  
C **  CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS  
C **  AT TIME LEVEL (N+1)  
C  
C----------------------------------------------------------------------C  
C 
      CALL CALCONC (ISTL,IS2TL)

C  
C----------------------------------------------------------------------C  
C  
      ! *** PMC BYPASS IF NOT SIMULATING SEDIMENTS
      IF(ISTRAN(6).GT.0.OR.ISTRAN(7).GT.0)THEN
        DO K=1,KB     
          DO L=1,LC  
            SEDBT(L,K)=0.  
            SNDBT(L,K)=0.  
          ENDDO  
        ENDDO  
C  
        DO NS=1,NSED  
          DO K=1,KB  
            DO L=1,LC  
              SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDDO  
C  
        DO NS=1,NSND  
          DO K=1,KB  
            DO L=1,LC  
              SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDDO  
C  
        DO K=1,KC  
          DO L=1,LC  
            SEDT(L,K)=0.  
            SNDT(L,K)=0.  
          ENDDO  
        ENDDO  
C  
        DO NS=1,NSED  
          DO K=1,KC  
            DO L=1,LC  
              SEDT(L,K)=SEDT(L,K)+SED(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDDO  
C  
        DO NS=1,NSND  
          DO K=1,KC  
            DO L=1,LC  
              SNDT(L,K)=SNDT(L,K)+SND(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF
C  
C----------------------------------------------------------------------C  
C  
C **  CHECK RANGE OF SALINITY AND DYE CONCENTRATION  
C  
      IF(ISMMC.EQ.1)THEN  
C  
        SALMAX=-100000.  
        SALMIN=100000.  
        DO K=1,KC  
          DO L=2,LA  
            IF(SAL(L,K).GT.SALMAX)THEN  
              SALMAX=SAL(L,K)  
              IMAX=IL(L)  
              JMAX=JL(L)  
              KMAX=K  
            ENDIF  
            IF(SAL(L,K).LT.SALMIN)THEN  
              SALMIN=SAL(L,K)  
              IMIN=IL(L)  
              JMIN=JL(L)  
              KMIN=K  
            ENDIF  
          ENDDO  
        ENDDO  
C  
        WRITE(6,6001)N  
        WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX  
        WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN  
C  
        SALMAX=-100000.  
        SALMIN=100000.  
        DO K=1,KC  
          DO L=2,LA  
            IF(DYE(L,K).GT.SALMAX)THEN  
              SALMAX=DYE(L,K)  
              IMAX=IL(L)  
              JMAX=JL(L)  
              KMAX=K  
            ENDIF  
            IF(DYE(L,K).LT.SALMIN)THEN  
              SALMIN=DYE(L,K)  
              IMIN=IL(L)  
              JMIN=JL(L)  
              KMIN=K  
            ENDIF  
          ENDDO  
        ENDDO  
C  
        WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX  
        WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN  
C  
        WRITE(8,6004)SALMAX,IMAX,JMAX,KMAX  
        WRITE(8,6005)SALMIN,IMIN,JMIN,KMIN  
C  
        SALMAX=-100000.  
        SALMIN=100000.  
        DO K=1,KC  
          DO L=2,LA  
            IF(SFL(L,K).GT.SALMAX)THEN  
              SALMAX=SFL(L,K)  
              IMAX=IL(L)  
              JMAX=JL(L)  
              KMAX=K  
            ENDIF  
            IF(SFL(L,K).LT.SALMIN)THEN  
              SALMIN=SFL(L,K)  
              IMIN=IL(L)  
              JMIN=JL(L)  
              KMIN=K  
            ENDIF  
          ENDDO  
        ENDDO  
C  
        WRITE(6,6006)SALMAX,IMAX,JMAX,KMAX  
        WRITE(6,6007)SALMIN,IMIN,JMIN,KMIN  
C  
      ENDIF  
C  
C  
      IF(ISMMC.EQ.2)THEN  
C  
        SALMAX=-100000.  
        SALMIN=100000.  
        DO K=1,KC  
          DO L=2,LA  
            IF(TEM(L,K).GT.SALMAX)THEN  
              SALMAX=TEM(L,K)  
              IMAX=IL(L)  
              JMAX=JL(L)  
              KMAX=K  
            ENDIF  
            IF(TEM(L,K).LT.SALMIN)THEN  
              SALMIN=TEM(L,K)  
              IMIN=IL(L)  
              JMIN=JL(L)  
              KMIN=K  
            ENDIF  
          ENDDO  
        ENDDO  
C  
        WRITE(6,6001)N  
        WRITE(6,6008)SALMAX,IMAX,JMAX,KMAX  
        WRITE(6,6009)SALMIN,IMIN,JMIN,KMIN  
C  
      ENDIF  
C  
 6001 FORMAT('  N=',I10)  
 6002 FORMAT('  SALMAX=',F14.4,5X,'I,J,K=',(3I10))  
 6003 FORMAT('  SALMIN=',F14.4,5X,'I,J,K=',(3I10))  
 6004 FORMAT('  DYEMAX=',F14.4,5X,'I,J,K=',(3I10))  
 6005 FORMAT('  DYEMIN=',F14.4,5X,'I,J,K=',(3I10))  
 6006 FORMAT('  SFLMAX=',F14.4,5X,'I,J,K=',(3I10))  
 6007 FORMAT('  SFLMIN=',F14.4,5X,'I,J,K=',(3I10))  
 6008 FORMAT('  TEMMAX=',F14.4,5X,'I,J,K=',(3I10))  
 6009 FORMAT('  TEMMIN=',F14.4,5X,'I,J,K=',(3I10))  

      ! *** DSLLC
      IF(DEBUG)THEN
        BTEST=.FALSE.
        LTEST=.FALSE.
        DO L=2,LA
          IF(ISNAN(HP(L)))THEN
            BTEST=.TRUE.
           IF(.NOT.LTEST)THEN
              OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',POSITION='APPEND',
     &             STATUS='UNKNOWN')  
              WRITE(1,*)' * DEBUG: ERROR IN DEPTH VARIABLES'
            ENDIF
            WRITE(1,910) TIMEDAY, L, IL(L), JL(L),
     &                   HP(L),H1P(L)
            HP(L)=H1P(L)
            LTEST=.TRUE.
          ENDIF
        ENDDO
        IF(LTEST)CLOSE(1,STATUS='KEEP')
        LTEST=.FALSE.
        IF(KC.GT.1)THEN
          DO L=2,LA
            DO K=1,KC
              IF(ISNAN(AV(L,K)))THEN
               BTEST=.TRUE.
                IF(.NOT.LTEST)THEN
                  OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',POSITION='APPEND',
     &                 STATUS='UNKNOWN')  
                  WRITE(1,*)' * DEBUG: ERROR IN VERTICAL VISCOSITY'
                ENDIF
                WRITE(1,9101) TIMEDAY, L, IL(L), JL(L), K, 'AV ',
     &                       AV(L,K)
                LTEST=.TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        IF(LTEST)CLOSE(1,STATUS='KEEP')

        LTEST=.FALSE.
        IF(ISTRAN(1).GE.1)THEN
          DO L=2,LA
            DO K=1,KC
              IF(ISNAN(SAL(L,K)))THEN
                BTEST=.TRUE.
               IF(.NOT.LTEST)THEN
                  OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',POSITION='APPEND',
     &                 STATUS='UNKNOWN')  
                  WRITE(1,*)' * DEBUG: ERROR IN SALINITY VARIABLES'
                ENDIF
                WRITE(1,911) TIMEDAY, L, IL(L), JL(L), K,
     &                       SAL(L,K),SAL1(L,K)  
                SAL(L,K)=SAL1(L,K)
                LTEST=.TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        IF(LTEST)CLOSE(1,STATUS='KEEP')

        LTEST=.FALSE.
        IF(ISTRAN(2).GE.1)THEN
          DO L=2,LA
            DO K=1,KC
              IF(ISNAN(TEM(L,K)))THEN
               BTEST=.TRUE.
                IF(.NOT.LTEST)THEN
                  OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',POSITION='APPEND',
     &                 STATUS='UNKNOWN')  
                  WRITE(1,*)' * DEBUG: ERROR IN TEMPERATURE VARIABLES'
                ENDIF
                WRITE(1,912) TIMEDAY, L, IL(L), JL(L), K,
     &                       TEM(L,K),TEM1(L,K)  
                TEM(L,K)=TEM1(L,K)
                LTEST=.TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        IF(LTEST)CLOSE(1,STATUS='KEEP')

        LTEST=.FALSE.
        IF(ISTRAN(6).GE.1)THEN
          DO L=2,LA
            DO K=1,KC
              DO NS=1,NSED
                IF(ISNAN(SED(L,K,NS)))THEN
                  BTEST=.TRUE.
                 IF(.NOT.LTEST)THEN
                    OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',POSITION='APPEND',
     &                 STATUS='UNKNOWN')  
                    WRITE(1,*)' * DEBUG: ERROR IN SED VARIABLES'
                  ENDIF
                  WRITE(1,916) TIMEDAY, L, IL(L), JL(L), K, NS, 
     &                         SED(L,K,NS),SED1(L,K,NS)  
                  SED(L,K,NS)=SED1(L,K,NS)
                  LTEST=.TRUE.
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF(LTEST)CLOSE(1,STATUS='KEEP')

        LTEST=.FALSE.
        IF(ISTRAN(7).GE.1)THEN
          DO L=2,LA
            DO K=1,KC
              DO NS=1,NSND
                IF(ISNAN(SND(L,K,NS)))THEN
                 BTEST=.TRUE.
                  IF(.NOT.LTEST)THEN
                    OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',POSITION='APPEND',
     &                   STATUS='UNKNOWN')  
                    WRITE(1,*)' * DEBUG: ERROR IN SAND VARIABLES'
                  ENDIF
                  WRITE(1,917) TIMEDAY, L, IL(L), JL(L), K, NS, 
     &                         SND(L,K,NS),SND(L,K,NS)  
                  SND(L,K,NS)=SND1(L,K,NS)
                  LTEST=.TRUE.
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF(LTEST)CLOSE(1,STATUS='KEEP')

        LTEST=.FALSE.
        IF(ISTRAN(8).GE.1)THEN
          DO L=2,LA
            DO K=1,KC
              DO NW=1,21
                IF(ISNAN(WQV(L,K,NW)))THEN
                 BTEST=.TRUE.
                  OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',POSITION='APPEND',
     &                 STATUS='UNKNOWN')  
                  WRITE(1,*)' * DEBUG: ERROR IN WATER QUALITY VARIABLES'
                  WRITE(1,918) TIMEDAY, L, IL(L), JL(L), K, NW, 
     &                         WQV(L,K,NW),WQVO(L,K,NW)  
                  CLOSE(1,STATUS='KEEP')
                  WQV(L,K,NW)=WQVO(L,K,NW)
                  LTEST=.TRUE.
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF(LTEST)CLOSE(1,STATUS='KEEP')
        
        ! *** DUMP THE RESULTS (JUST PRIOR) TO EE FOR ANALYSIS
        IF(BTEST)THEN
          CALL SURFPLT  
          CALL VELPLTH
          CALL EEXPOUT(-1)  
          CLOSE(7)  
          CLOSE(8)  
          CLOSE(9)  
          STOP 'ERROR: NANs have been computed!'
        ENDIF
  910 FORMAT('ERROR: TIME, L, I, J, HP = ',        F10.5,3I6,2F10.4)
 9101 FORMAT('ERROR: TIME, L, I, J, K, ',A3,' = ', F10.5,4I6,2F10.4)
  911 FORMAT('ERROR: TIME, L, I, J, K, SAL = ',    F10.5,4I6,2F10.4)
  912 FORMAT('ERROR: TIME, L, I, J, K, TEM = ',    F10.5,4I6,2F10.4)
  916 FORMAT('ERROR: TIME, L, I, J, K, NS, SED = ',F10.5,5I6,2F10.4)
  917 FORMAT('ERROR: TIME, L, I, J, K, NX, SND = ',F10.5,5I6,2F10.4)
  918 FORMAT('ERROR: TIME, L, I, J, K, NW, WQV = ',F10.5,5I6,2F10.4)
      ENDIF
C  
! DEBUG NAN IN DYE FEARGHAL 
      BTEST = .FALSE.
      DO K = 1,KC
         DO L =2,LA
          IF(ISNAN(DYE(L,K)))THEN
                  write(*,*) 'DYE(L,K)',DYE(L,K),XPAR(IL(L)),YPAR(JL(L))
                  BTEST=.TRUE.
             I = XPAR(IL(L)); J = YPAR(JL(L))
            OPEN(1,FILE='ERROR'//ans(partid2)//'.LOG',status='unknown',position='append')
            write(1,*) I,J,K,N,LMASKDRY(L),IMASKDRY(L),DYE(L,K),HP(L),U(L,K),V(L,K),
     &       SUB(L),SUB(LEAST(L)),SVB(L),SVB(LN),HP(LEAST(L)),UHDYE(LEAST(L))   
            close(1)
          END IF
        END DO
      END DO
       IF(BTEST)THEN
          CALL SURFPLT  
          CALL VELPLTH
          CALL TMSR  
          CLOSE(7)  
          CLOSE(8)  
          CLOSE(9)  
          STOP 'ERROR: NANs have been computed!'
        ENDIF


C**********************************************************************C  
C  
C **  CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT  
C **  CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOUBLE TIME  
C **  STEP TRANSPORT FIELD  
C  
C----------------------------------------------------------------------C  
C  
      ITMP=0  
      IF(ISTRAN(4).GE.1) ITMP=1  
      IF(ISTRAN(8).GE.1) ITMP=1  
      IF(ISWASP.GE.1)ITMP=1              ! 6/7/2005 a stoddard dsllc 
      IF(ISICM.GE.1) ITMP=1  
C  
      IF(ITMP.EQ.1)THEN  
C  
C **  CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION  
C  
        DO L=2,LA  
          HWQ(L)=HP(L)  
          WWQ(L,0)=0.  
        ENDDO  
C  
        DO K=1,KC  
          DO L=2,LA  
            UHDYWQ(L,K)=UHDY2(L,K)  
            VHDXWQ(L,K)=VHDX2(L,K)  
            UWQ(L,K)=U2(L,K)  
            VWQ(L,K)=V2(L,K)  
            WWQ(L,K)=W2(L,K)  
          ENDDO  
        ENDDO  
C  
C     ADD CHANNEL INTERACTIONS  
C  
  
        IF(MDCHH.GE.1)THEN  
          DO NMD=1,MDCHH  
            IF(MDCHTYP(NMD).EQ.1)THEN  
              HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))  
     &            +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))  
              HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))  
     &            -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))  
            ENDIF  
            IF(MDCHTYP(NMD).EQ.2)THEN  
              HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))  
     &            +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))  
              HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))  
     &            -DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))  
            ENDIF  
            IF(MDCHTYP(NMD).EQ.3)THEN  
              HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))  
     &            +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))  
     &            +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))  
              HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))  
     &            -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))  
              HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))  
     &            -DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))  
            ENDIF  
          ENDDO  
        ENDIF  
C  
C     END ADD CHANNEL INTERACTIONS  
C  
        IF(ISTRAN(8).GE.1) CALL WQ3D(ISTL,IS2TL)  
        IF(ISTRAN(4).GE.1) CALL CALSFT(ISTL,IS2TL)
C  
        DO L=2,LA  
          H2WQ(L)=HWQ(L)  
        ENDDO  
C  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING  
C **  AN EQUATION OF STATE  
C  
      DO K=1,KC  
        DO L=2,LA  
          B1(L,K)=B(L,K)  
        ENDDO  
      ENDDO  
C  
      IF(BSC.GT.1.E-6)THEN  
        CALL CALBUOY  
      ELSE  
        DO K=1,KC  
          DO L=2,LA  
            B(L,K)=0.  
          ENDDO  
        ENDDO  
      ENDIF  
C  
C **  CALL TWO-TIME LEVEL BALANCES  
C  
      IF(ISBAL.GE.1)THEN  
        CALL BAL2T4  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)  
C  
C----------------------------------------------------------------------C  
C  
#ifdef key_mpi
        CALL COMMUNICATE_P(U(:,1))
        CALL COMMUNICATE_P(V(:,1))
#endif
      DO L=2,LA  
        LN=LNC(L)  
        LS=LSC(L)  
        LE=LEAST(L)
        LW=LWEST(L)
        LNW=LNWC(L)  
        LSE=LSEC(L)  
        LSW=LSWC(L)  
        U1V(L)=UV(L)  
        V1U(L)=VU(L)  
        UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))  
     &             +HP(L )*(U(LE ,1)+U(L ,1)))*HVI(L)  
        VU(L)=0.25*(HP(LW)*(V(LNW,1)+V(L,1))  
     &             +HP(L )*(V(LN ,1)+V(L,1)))*HUI(L)  
      ENDDO  
C  
C**********************************************************************C  
C  
C **  CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES  
C **  AT TIME LEVEL (N)  
C  
      IF(ISHDMF.GE.1) CALL CALHDMF
      
C  
C**********************************************************************C  
C  
C **  CALCULATE BOTTOM STRESS AT LEVEL (N+1)  
C  
      T1TMP=SECNDS(SECNDS_ZERO)  
C
      CALL CALTBXY(ISTL,IS2TL)  
C  
      DO L=2,LA  
        TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)  
     &      +U(L,1)*U(L,1)))*U(L,1)  
        TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)  
     &      +V(L,1)*V(L,1)))*V(L,1)  
      ENDDO  
C  
C**********************************************************************C  
C  
C **  SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES  
C  
      IF(ISBSDFUF.GE.1)THEN  
        HDFUFM=1.E-12  
C  
        DO L=2,LA  
          LS=LSC(L)  
          HDFUFX(L)=HDFUFM+G*SUB(L)*HU(L)*(BELV(LWEST(L))-BELV(L))*DXIU(L)  
          HDFUFY(L)=HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )-BELV(L))*DYIV(L)  
        ENDDO  
C  
        DO L=2,LA  
          HDFUFX(L)=TBX(L)/HDFUFX(L)  
          HDFUFY(L)=TBY(L)/HDFUFY(L)  
        ENDDO  
C  
        DO L=2,LA  
          HDFUFX(L)=MAX(HDFUFX(L),-1.0)  
          HDFUFY(L)=MAX(HDFUFY(L),-1.0)  
        ENDDO  
C  
        DO L=2,LA  
          HDFUFX(L)=MIN(HDFUFX(L),1.0)  
          HDFUFY(L)=MIN(HDFUFY(L),1.0)  
        ENDDO  
C  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)  
C  
C----------------------------------------------------------------------C  
C  
      IF(ISWAVE.EQ.0)THEN  
C  
C----------------------------------------------------------------------c  
C  
        IF(ISCORTBC.EQ.0) THEN  
C  
          DO L=2,LA  
            TVAR3S(L)=TSY(LNC(L))  
            TVAR3W(L)=TSX(LEAST(L))  
            TVAR3E(L)=TBX(LEAST(L)   )  
            TVAR3N(L)=TBY(LNC(L))  
          ENDDO  
C  
          DO L=2 ,LA  
            TMP = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2  
     &           +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2  
            QQ(L,0 )=0.5*CTURB2*SQRT(TMP)

            TMP = (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2  
     &           +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2  
            QQ(L,KC)=0.5*CTURB2*SQRT(TMP) 
             
            QQSQR(L,0)=SQRT(QQ(L,0))  ! *** DSLLC
          ENDDO  
C  
        ENDIF  
C  
C----------------------------------------------------------------------c  
C  
        IF(ISCORTBC.GE.1) THEN  
C  
          IF(ISCORTBCD.GE.1)THEN  
            NTMPVAL=MOD(N,NTSPTC)  
            IF(NTMPVAL.EQ.0.AND.DEBUG)THEN  
              OPEN(1,FILE='ADJSTRESSE.OUT',POSITION='APPEND')  
            ENDIF  
          ENDIF  
C  
          DO L=2,LA  
            TVAR3S(L)=TSY(LNC(L))  
            TVAR3W(L)=TSX(LEAST(L))  
            TVAR3E(L)=TBX(LEAST(L)   )  
            TVAR3N(L)=TBY(LNC(L))  
          ENDDO  
C  
          DO L=2,LA  
            WCOREST(L)=1.  
            WCORWST(L)=1.  
            WCORNTH(L)=1.  
            WCORSTH(L)=1.  
          ENDDO  
C  
          DO L=2,LA  
            IF(ISSBCP(L).EQ.0)THEN  
              IF(SUB(LEAST(L)).LT.0.5)WCOREST(L)=FSCORTBCV(L)  
              IF(SUB(L).LT.0.5)WCORWST(L)=FSCORTBCV(L)  
              IF(SVB(LNC(L)).LT.0.5)WCORNTH(L)=FSCORTBCV(L)  
              IF(SVB(L).LT.0.5)WCORSTH(L)=FSCORTBCV(L)  
            ENDIF  
          ENDDO  
C  
          DO L=2,LA  
            WCOREW(L)=1./(WCOREST(L)+WCORWST(L))  
            WCORNS(L)=1./(WCORNTH(L)+WCORSTH(L))  
          ENDDO  
C  
          DO L=2,LA  
            WCOREST(L)=WCOREST(L)*WCOREW(L)  
            WCORWST(L)=WCORWST(L)*WCOREW(L)  
            WCORNTH(L)=WCORNTH(L)*WCORNS(L)  
            WCORSTH(L)=WCORSTH(L)*WCORNS(L)  
          ENDDO  
C  
          DO L=2,LA  
            QQ(L,0 )=CTURB2*SQRT(  
     &          (RSSBCE(L)*WCOREST(L)*TVAR3E(L)  
     &          +RSSBCW(L)*WCORWST(L)*TBX(L))**2  
     &          +(RSSBCN(L)*WCORNTH(L)*TVAR3N(L)  
     &          +RSSBCS(L)*WCORSTH(L)*TBY(L))**2)  
            QQ(L,KC)=0.5*CTURB2*SQRT(  
     &          (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2  
     &          +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)  
            QQSQR(L,0)=SQRT(QQ(L,0))  ! *** DSLLC
          ENDDO  
C  
          IF(ISCORTBCD.GE.1.AND.NTMPVAL.EQ.0)THEN  
C  
            DO L=2,LA  
              LCORNER(L)=0  
              KCORNER=0  
              IF(WCORWST(L).GT.0.505)THEN  
                KCORNER=KCORNER+1  
                LCORNWE(L)=LWEST(L)  
              ENDIF  
              IF(WCOREST(L).GT.0.505)THEN  
                KCORNER=KCORNER+1  
                LCORNWE(L)=LEAST(L)  
              ENDIF  
              IF(WCORNTH(L).GT.0.505)THEN  
                KCORNER=KCORNER+1  
                LCORNSN(L)=LNC(L)  
              ENDIF  
              IF(WCORSTH(L).GT.0.505)THEN  
                KCORNER=KCORNER+1  
                LCORNSN(L)=LSC(L)  
              ENDIF  
              IF(KCORNER.EQ.2)LCORNER(L)=1  
            ENDDO  
C  
            NCORCELLS=0  
            DO L=2,LA  
              NCORCELLS=NCORCELLS+LCORNER(L)  
            ENDDO  
C  
            IF(DEBUG)THEN
              WRITE(1,3675)TIMEDAY,NCORCELLS  
C  
              DO L=2,LA  
                IF(LMASKDRY(L))THEN  
                  IF(LCORNER(L).EQ.1)THEN  
                    LWE=LCORNWE(L)  
                    LSN=LCORNSN(L)  
                    TAUTMP=QQ(L,0)/CTURB2  
                    TAUTMPWE=QQ(LWE,0)/CTURB2  
                    TAUTMPSN=QQ(LSN,0)/CTURB2  
                    WRITE(1,3677)IL(L),JL(L),TAUTMP,TAUBSND(L),  
     &                TAUBSED(L)  
                    WRITE(1,3676)IL(LWE),JL(LWE),TAUTMPWE,TAUBSND(LWE),  
     &                TAUBSED(LWE)  
                    WRITE(1,3676)IL(LSN),JL(LSN),TAUTMPSN,TAUBSND(LSN),  
     &                TAUBSED(LSN)  
                  ENDIF  
                ENDIF  
              ENDDO  
            ENDIF
C  
          ENDIF  
  
          IF(DEBUG)CLOSE(1)  
C  
        ENDIF  
C  
C----------------------------------------------------------------------c  
C  
      ENDIF  
C  
 3678 FORMAT(2I6,4F13.3)  
 3679 FORMAT(12x,4F13.3)  
 3680 FORMAT(12x,6F13.5)  
 3681 FORMAT(12X,5E13.4,F13.5)  
 3677 FORMAT('CORNER',2I5,5E14.5)  
 3676 FORMAT(6X,2I5,5E14.5)  
 3675 FORMAT(F11.3,I6,' TIME IN DAYS AND NUMBER OF CORNERS')  
C  
C     ENDIF  
C  
C**********************************************************************C  
C  
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)  
C  
C----------------------------------------------------------------------C  
C  
      IF(ISWAVE.GE.1)THEN  
C  
        DO L=2,LA  
          TVAR3S(L)=TSY(LNC(L))  
          TVAR3W(L)=TSX(LEAST(L))  
          TVAR3E(L)=TBX(LEAST(L)   )  
          TVAR3N(L)=TBY(LNC(L))  
        ENDDO  
C  
        DO L=2,LA  
          TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2  
     &            +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2  
          TAUBC=0.5*SQRT(TAUBC2)  
          UTMP=0.5*STCUV(L)*(U(LEAST(L),1)+U(L,1))+1.E-12  
          VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))  
          CURANG=ATAN2(VTMP,UTMP)  
          TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))  
     &        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))  
          TAUB2=MAX(TAUB2,0.)  
          QQ(L,0 )=CTURB2*SQRT(TAUB2)  
          QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2  
     &        +(TVAR3S(L)+TSY(L))**2)  
          QQSQR(L,0)=SQRT(QQ(L,0))  ! *** DSLLC
        ENDDO  
C  
      ENDIF  
C
      TTBXY=TTBXY+SECNDS(T1TMP)  
C  
C**********************************************************************C  
C  
C **  CALCULATE TURBULENT INTENSITY SQUARED  
C   
      T1TMP=SECNDS(SECNDS_ZERO)  
      IF(KC.GT.1)THEN  
        IF(ISQQ.EQ.1)THEN  
          IF(ISTOPT(0).EQ.0)CALL CALQQ2TOLD (ISTL)  
          IF(ISTOPT(0).GE.1)CALL CALQQ2T (ISTL) 
        ELSEIF(ISQQ.EQ.2)THEN
          CALL CALQQ2 (ISTL)  
        ENDIF
      ENDIF  
      TQQQ=TQQQ+SECNDS(T1TMP)  
C  
C**********************************************************************C  
C  
C **  CALCULATE MEAN MASS TRANSPORT FIELD  
C  
      IF(ISSSMMT.NE.2)THEN  
        IF(ISICM.GE.1)THEN  
          NTMP=MOD(N,2)  
          IF(ISTL.EQ.3.AND.NTMP.EQ.0) CALL CALMMT  
        ENDIF  
      ENDIF  
C  
C      IF(ISSSMMT.NE.2) CALL CALMMT  
C  
C**********************************************************************C  
C  
C **  HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED  
C  
C**********************************************************************C  
C**********************************************************************C  
C  
C **  UPDATE SOLUTION WITH DATA ASSIMILATION INNOVATION  
C  
C**********************************************************************C  
      DAFLAG = 0   ! include flag in input files for data assimilation 0/1
      IF (DAFLAG == 1) THEN ! make DA call  
        IF (TIMESEC .GE. TSTARTSEC .AND.   ! Time flags must be
     &      TIMESEC .LE. TENDSEC) THEN     ! correctly implemented
!          CALL BLUE_INI
         END IF   
      END IF  

C**********************************************************************C
C
C **  CALL FOR DATA ASSIMILATION OF VERTICAL PROFILER DATA
C
C**********************************************************************C
#ifdef key_da
      IF (IDA_FLAG == 1) THEN ! DATA ASSIMILATION ACTIVE
        IF (N == NDATA_ASSIM) THEN ! IF TIME = DATA ASSIMILATION STEP BASED ON FREQ
          CALL DA_INI
          NDATA_ASSIM = NDATA_ASSIM + (DAFREQ * 3600/DT)
         END IF
      END IF
#endif
C**********************************************************************C
C
C     communicate variables with MPI calls before writing to files and looping computations
C
C     
#ifdef key_mpi
        CALL COMMUNICATE_MPI 
        IF(ISTRAN(1).EQ.1)  CALL COMMUNICATE_3D(SAL)
        IF(ISTRAN(2).EQ.1)  CALL COMMUNICATE_3D(TEM)
        IF(ISTRAN(3).EQ.1)  CALL COMMUNICATE_3D(DYE)
        IF(ISTRAN(4).EQ.1)  CALL COMMUNICATE_3D(SFL)
        IF(ISTRAN(5).EQ.1)  CALL COMMUNICATE_3D(TVAR1S)
        IF(ISTRAN(6).EQ.1)  CALL COMMUNICATE_3D(SEDT)
        IF(ISTRAN(7).EQ.1)  CALL COMMUNICATE_3D(SNDT)
        IF(ISTRAN(8).EQ.1)  CALL COMMUNICATE_4D_WQ(WQV)
#endif
C**********************************************************************C
C
C **  OPTIONAL WE DO A CHECK ON MASS CONSERVATION FOR FRESHWATER INFLOWS
C
C**********************************************************************C
      IF (KINVARA_MC) THEN
          TOTVOL=0.0
          TOTVOL0=0.0
          TOTVOL00=0.0
          VFW0=0.0
          VSW=0.0
          VFW00=0.0
          VFW1=0.0
          DO I = 3, IC-2
            DO J = 3, JC - 2
            IF (IJCT(I,J) == 5) THEN
               L = LIJ(I,J)
                TOTVOL=TOTVOL+((.5*(HU(L) + HV(L)))*DXYP(L))
                SAL_CELL = SUM(SAL(L,:))/(1. * KC)
                VFW1= VFW1+((.5*( HU(L) + HV(L)))*DXYP(L)*(1-( SAL_CELL/35)))
                VSW= VSW+((.5*( HU(L) + HV(L)))*DXYP(L)*(SAL_CELL/35))
             END IF
             END DO
          END DO
          VFW2= VFW2  + (SUM(QSUME)*DT)
!            write(*,*) QSS(1,1), DXP(10), HU(10), HV(10), BELV(10), IL(10), JL(10), DT
!     c              write(*,*) totvol0,totvol00,vsw0,vsw00,vfw0,vfw00
          IFLPRINT = IFLPRINT + 1
#ifdef key_mpi
          CALL MPI_ALLREDUCE(TOTVOL,TOTVOL_OUT,1,MPI_REAL8,MPI_SUM,EFDC_COMM,IERR)
          CALL MPI_ALLREDUCE(VFW1,VFW1_OUT,1,MPI_REAL8,MPI_SUM,EFDC_COMM,IERR)
          CALL MPI_ALLREDUCE(VFW2,VFW2_OUT,1,MPI_REAL8,MPI_SUM,EFDC_COMM,IERR)
          CALL MPI_ALLREDUCE(VSW,VSW_OUT,1,MPI_REAL8,MPI_SUM,EFDC_COMM,IERR)
#else
          TOTVOL_OUT = TOTVOL
          VFW1_OUT = VFW1
          VFW2_OUT = VFW2
          VSW_OUT = VSW
#endif


          if (IFLPRINT.eq.50) then
                write(777,777) timesec/86400,TOTVOL_OUT,(TOTVOL_OUT-VFW1_OUT-VFW2_OUT),VFW1_OUT, VSW_OUT,VFW2_OUT, VSW_BOUNDARY
                IFLPRINT = 0
          end if
      END IF
!      L = LIJ(XLOC(248), XLOC(29))
!      SAL_CELL = SUM(SAL(L,:))/KC
!      VFW_TEMP= ((1.*( HP(L)))*DXYP(L)*(1-( SAL_CELL/35)))
!      write(*,*) XLOC(248), XLOC(29), VFW_TEMP, HU(L), HV(L), HP(L)

 777  FORMAT(F13.5, 7E13.4)
C**************************************************************

C  
C **  WRITE TO TIME SERIES FILES  
C  
      IF(ISDYNSTP.EQ.0)THEN  
        CTIM=DT*FLOAT(N)+TCON*TBEGIN  
        CTIM=CTIM/TCON  
      ELSE  
        CTIM=TIMESEC/TCON  
      ENDIF  
C  
CDYN      IF(ISTMSR.GE.1)THEN  
CDYN        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN  
CDYN          IF(NCTMSR.EQ.NWTMSR)THEN  
CDYN            CALL TMSR  
CDYN            ICALLTP=1  
CDYN            NCTMSR=1  
CDYN           ELSE  
CDYN            NCTMSR=NCTMSR+1  
CDYN          ENDIF  
CDYN        ENDIF  
CDYN      ENDIF  
C  
C  
      IF(ISTMSR.GE.1)THEN  
c        IF(N.GE.NBTMSR.AND.N.LE.NSTMSR)THEN  
        IF(NCTMSR.GE.NWTMSR)THEN  
          CALL TMSR  
          NDIFF=NWTMSR-NCTMSR  
          ICALLTP=1  
          NCTMSR=NINCRMT+NDIFF  
        ELSE  
          NCTMSR=NCTMSR+NINCRMT  
        ENDIF  
c        ENDIF  
      ENDIF  
C  
C****************************************************
C **  WRITE TO DUMP FILES  ******************C  
C  
C  
      IF(ISDUMP.GE.1)THEN  
        IF(CTIM.GE.TSDUMP.AND.CTIM.LE.TEDUMP)THEN  
C          IF(NCDUMP.EQ.NSDUMP)THEN  
          IF(NCDUMP.GE.NSDUMP)THEN  
            CALL DUMP  
            NDIFF=NSDUMP-NCDUMP  
            ICALLTP=1  
C            NCDUMP=1  
            NCDUMP=NINCRMT+NDIFF  
          ELSE  
C            NCDUMP=NCDUMP+1  
            NCDUMP=NCDUMP+NINCRMT  
          ENDIF  
        ENDIF  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  OUTPUT ZERO DIMENSION VOLUME BALANCE  
C  
C----------------------------------------------------------------------C  
C  
      IF(ISDRY.GE.1.AND.ISDRY.LT.98)THEN  
        IF(ICALLTP.EQ.1.AND.DEBUG)THEN  
          OPEN(1,FILE='ZVOLBAL.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          DO LS=1,LORMAX  
            IF(VOLZERD.GE.VOLSEL(LS).AND.VOLZERD.LT.VOLSEL(LS+1))THEN  
              WTM=VOLSEL(LS+1)-VOLZERD  
              WTMP=VOLZERD-VOLSEL(LS)  
              DELVOL=VOLSEL(LS+1)-VOLSEL(LS)  
              WTM=WTM/DELVOL  
              WTMP=WTMP/DELVOL  
              SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)  
              ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)  
            ENDIF  
          ENDDO  
          IF(ISDYNSTP.EQ.0)THEN  
            CTIM=DT*FLOAT(N)+TCON*TBEGIN  
            CTIM=CTIM/TCTMSR  
          ELSE  
            CTIM=TIMESEC/TCTMSR  
          ENDIF  
          WRITE(1,5304) CTIM,SELZERD,ASFZERD,VOLZERD,VETZERD  
          CLOSE(1)  
        ENDIF  
      ENDIF  
      ICALLTP=0  
C  
 5304 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))  
C  
C**********************************************************************C  
C  
C **  WRITE VERTICAL SCALAR FIELD PROFILES  
C  
      IF(ISVSFP.EQ.1)THEN  
        IF(N.GE.NBVSFP.AND.N.LE.NSVSFP)THEN  
          CALL VSFP  
        ENDIF  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  CALCULATE MEAN MASS TRANSPORT FIELD  
C  
      IF(ISSSMMT.NE.2)THEN  
        IF(ISICM.EQ.0) CALL CALMMT  
      ENDIF  
C  
C      IF(ISSSMMT.NE.2) CALL CALMMT  
C  
C**********************************************************************C  
C  
C **  ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES  
C  
      IF(ISPD.GE.2)THEN                     !DHC
        IF (TIMEDAY.GE.LA_BEGTI.AND.TIMEDAY.LE.LA_ENDTI) THEN
          T1TMP=SECNDS(SECNDS_ZERO)                
          CALL DRIFTERC
          TLRPD=TLRPD+SECNDS(T1TMP)  
        ENDIF
      ENDIF  
      
!      IF(ISLRPD.GE.1)THEN  
!        T1TMP=SECNDS(SECNDS_ZERO)                  !DHC:13-04-09
!        IF(ISLRPD.LE.2)THEN  
!          IF(N.GE.NLRPDRT(1)) CALL LAGRES  
!        ENDIF  
!        IF(ISLRPD.GE.3)THEN  
!          IF(N.GE.NLRPDRT(1)) CALL GLMRES  
!        ENDIF  
!        TLRPD=TLRPD+SECNDS(T1TMP)  
!      ENDIF  
C  
C**********************************************************************C  
C  
C **  CALCULATE VOLUME MASS, MOMENTUM AND ENERGY BALANCES  
C  
C      IF(ISBAL.GE.1)THEN  
C         CALL CALBAL5  
C         NTMP=MOD(N,2)  
C         IF(NTMP.EQ.0)THEN  
C           CALL CBALEV5  
C          ELSE  
C           CALL CBALOD5  
C         ENDIF  
C       ENDIF  
C  
C   SEDIMENT BUDGET CALCULATION     (DLK 10/15)  
C  
C       IF(ISSBAL.GE.1)THEN  
C       CALL BUDGET5  
C       ENDIF  
C       NTMP=MOD(N,2)  
C       IF(NTMP.EQ.0)THEN  
C         CALL BUDGEV5  
C        ELSE  
C         CALL BUDGOD5  
C       ENDIF  
C  
C **  CALL TWO-TIME LEVEL BALANCES  
C  
      IF(ISBAL.GE.1)THEN  
        CALL BAL2T5  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  PERFORM AN M2 TIDE HARMONIC ANALYSIS EVERY 2 M2 PERIODS  
C  
      IF(ISHTA.EQ.1) CALL CALHTA  
C  
C**********************************************************************C  
C  
C **  CALCULATE DISPERSION COEFFICIENTS  
C  
C     IF(N.GE.NDISP)THEN  
      IF(N.GE.NDISP.AND.NCTBC.EQ.1)THEN  
        IF(ISDISP.EQ.2) CALL CALDISP2  
        IF(ISDISP.EQ.3) CALL CALDISP3  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  PERFORM LEAST SQUARES HARMONIC ANALYSIS AT SELECTED LOCATIONS  
C  
      IF(ISLSHA.EQ.1.AND.N.EQ.NCLSHA)THEN  
        CALL LSQHARM  
        NCLSHA=NCLSHA+(NTSPTC/24)  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  PRINT INTERMEDIATE RESULTS  
C  
C----------------------------------------------------------------------C  
C  
      IF(NPRINT .EQ. NTSPP)THEN  
        NPRINT=1  
        CALL OUTPUT1  
      ELSE  
        NPRINT=NPRINT+1  
      ENDIF  

C**********************************************************************C  
C  
C **  WRITE TO TIME VARYING GRAPHICS FILES  
C  
C----------------------------------------------------------------------C  
C  
      IPLTTMP=0
      IF(ISPPH.EQ.1.OR.ISPPH.EQ.2)IPLTTMP=1
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
Cpmc      IF(N.GE.NCPPH.AND.ISPPH.GE.1)THEN  
          CALL SURFPLT  
          NCPPH = NCPPH + (NTSPTC/NPPPH)
      END IF  
C  
C----------------------------------------------------------------------C  
C  
CDYN      IF(N.EQ.NCBPH.AND.ISBPH.EQ.1)THEN  
C  
      IF(N.GE.NCBPH.AND.ISBPH.GE.1)THEN  
        IF(ISBEXP.EQ.0)THEN  
          CALL BEDPLTH  
          NCBPH=NCBPH+(NTSPTC/NPBPH)  
        ENDIF  
      ENDIF  
C  
C----------------------------------------------------------------------C  
C  
CDYN      IF(N.EQ.NCVPH.AND.ISVPH.GE.1)THEN  
C        IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN  

      IPLTTMP=0  
      IF(ISVPH.EQ.1.OR.ISVPH.EQ.2)IPLTTMP=1 
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
!      IF(N.EQ.NCVPH.AND.IPLTTMP==1)THEN
        CALL VELPLTH 
!        NCVPH=NCVPH+(NTSPTC/NPVPH)
      ENDIF  

C  
C----------------------------------------------------------------------C  
C  
CDYN      IF(N.EQ.NCVPV.AND.ISVPV.GE.1)THEN  
C  
      IF(N.GE.NCVPV.AND.ISVPV.GE.1)THEN  
        CALL VELPLTV  
        NCVPV=NCVPV+(NTSPTC/NPVPV)  
      ENDIF  
C  
C----------------------------------------------------------------------C  
C  
      DO K=1,KC  
        DO L=1,LC  
          TVAR1S(L,K)=TOX(L,K,1)  
        ENDDO  
      ENDDO  
C  

      IPLTTMP=0  
      IF(ISSPH(1).EQ.1.OR.ISSPH(1).EQ.2)IPLTTMP=1  
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
        IF(ISTRAN(1).GE.1) CALL SALPLTH (1,SAL)  
        NCSPH(1)=NCSPH(1)+(NTSPTC/NPSPH(1))  
      ENDIF  
C  
      IPLTTMP=0  
      IF(ISSPH(2).EQ.1.OR.ISSPH(2).EQ.2)IPLTTMP=1  
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
        IF(ISTRAN(2).GE.1) CALL SALPLTH (2,TEM)  
        NCSPH(2)=NCSPH(2)+(NTSPTC/NPSPH(2))  
      ENDIF  
C  
      IPLTTMP=0  
      IF(ISSPH(3).EQ.1.OR.ISSPH(3).EQ.2)IPLTTMP=1  
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
        IF(ISTRAN(3).GE.1) CALL SALPLTH (3,DYE)  
        NCSPH(3)=NCSPH(3)+(NTSPTC/NPSPH(3))  
      ENDIF  
C  
      IPLTTMP=0  
      IF(ISSPH(4).EQ.1.OR.ISSPH(4).EQ.2)IPLTTMP=1  
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
        IF(ISTRAN(4).GE.1) CALL SALPLTH (4,SFL)  
        NCSPH(4)=NCSPH(4)+(NTSPTC/NPSPH(4))  
      ENDIF  
C  
      IPLTTMP=0  
      IF(ISSPH(5).EQ.1.OR.ISSPH(5).EQ.2)IPLTTMP=1  
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
        IF(ISTRAN(5).GE.1) CALL SALPLTH (5,TVAR1S)  
        NCSPH(5)=NCSPH(5)+(NTSPTC/NPSPH(5))  
      ENDIF  
C  
      IPLTTMP=0  
      IF(ISSPH(6).EQ.1.OR.ISSPH(6).EQ.2)IPLTTMP=1  
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
        IF(ISTRAN(6).GE.1) CALL SALPLTH (6,SEDT)  
        NCSPH(6)=NCSPH(6)+(NTSPTC/NPSPH(6))  
      ENDIF  
C  
      IPLTTMP=0  
      IF(ISSPH(7).EQ.1.OR.ISSPH(7).EQ.2)IPLTTMP=1  
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS).AND.IPLTTMP.EQ.1)THEN
        IF(ISTRAN(7).GE.1) CALL SALPLTH (7,SNDT)  
        NCSPH(7)=NCSPH(7)+(NTSPTC/NPSPH(7))  
      ENDIF  
C  
C----------------------------------------------------------------------C  
C  
      DO ITMP=1,7  
        IF(N.GE.NCSPV(ITMP).AND.ISSPV(ITMP).GE.1)THEN  
          CALL SALPLTV(ITMP)  
          NCSPV(ITMP)=NCSPV(ITMP)+(NTSPTC/NPSPV(ITMP))  
        ENDIF  
      ENDDO  
C  
C----------------------------------------------------------------------C  
C  
C **  WRITE EFDC EXPLORER FORMAT OUTPUT  
C  
      IF(ISSPH(8).EQ.1.AND.ISBEXP.EQ.1)THEN  
        IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS))THEN  
          CALL EEXPOUT(0)  
        ENDIF  
      ENDIF  
      IF(TIMEDAY.GE.SNAPSHOTS(NSNAPSHOTS))THEN  
        NSNAPSHOTS=NSNAPSHOTS+1
      ENDIF  
C  
C**********************************************************************C  
C  
C **  WRITE TO TIME VARYING 3D HDF GRAPHICS FILES  
C  
C----------------------------------------------------------------------C  
C  
      IF(N.EQ.NC3DO.AND.IS3DO.EQ.1)THEN  
        CALL OUT3D  
        NC3DO=NC3DO+(NTSPTC/NP3DO)  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  WRITE RESTART FILE EVERY ISRESTO M2 TIDAL CYCLES  
C  
      IF(ISRESTO.GE.1)THEN  
        IF((N-ISSREST).GT.NRESTO)THEN  
          CALL RESTOUT(0)  
          IF(ISTRAN(8).GE.1)THEN  
            IF(IWQRST.EQ.1) CALL WWQRST  
            IF(IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST  
          ENDIF
          ISSREST=N  
        ENDIF  
      ENDIF
      
C  
C**********************************************************************C  
C  
C **  RECORD TIME  
C  
C **  DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE  
C **  SUPPORTED ON OTHER SYSTEMS.  
C  
      IF(NTIMER.EQ.NTSPTC)THEN  
C *** EE BEGIN BLOCK  
        IF (PARTID == 0 ) CALL TIMELOG(N,TIMEDAY)
C *** EE END BLOCK  
        NTIMER=1  
      ELSE  
        NTIMER=NTIMER+1  
      ENDIF  

C  
C**********************************************************************C  
C  

      IF(ISHOW.GT.0)CALL SHOWVAL

C  
C**********************************************************************C  
C  
C *** DJB  
!      IF(.NOT.KBHIT())GOTO 1001  
!      I1=GETCH()  
!      WRITE(*,'(A)')'PROGRAM PAUSED BY USER'  
!      WRITE(*,'(A)')'  EFDC_DS: TO EXIT PRESS THE SAME KEY'
!      WRITE(*,'(A)')'  EFDC_DS: TO CONTINUE RUN PRESS ANY OTHER KEY'  
!      I2=GETCH()  
 !     IF(I1.NE.I2)GOTO 1001
C  
      DO L = 2,LA
        PPTMP=HP(L)+BELV(L)   ! tidal elevation from MWL
        MXR(L) = MAX(MXR(L),PPTMP)
        MNR(L) = MIN(MNR(L),PPTMP)
      END DO


      GOTO 1001
 1000 CONTINUE  

      IF (DEBUG) THEN
        IF (PARTID == MASTER_TASK) THEN
          OPEN(112,FILE='TIMING.DAT',STATUS='UNKNOWN',POSITION='APPEND')
          CALL CPU_TIME(OEND)
          WRITE(112,*) NPARTS,OEND-OSTART
          CLOSE(112)
        END IF
      ENDIF
C  
C**********************************************************************C  
C  
C **  TIME LOOP COMPLETED  
C  
      THDMT=THDMT+SECNDS(TTMP)  
C  
C**********************************************************************C  
C *** EE BEGIN BLOCK  
C       MOVED THE TIMING OUTPUT BLOCK TO THE MAIN AAEFDC TO ELIMINATE  
C       UNNECESSARY DUPLICATION  
C *** EE END BLOCK  
C**********************************************************************C  
C  
 2000 CONTINUE  
C  
C**********************************************************************C  
C  
C **  PRINT FINAL RESULTS  
C  
      CALL OUTPUT2  

C  
**********************************************************************C  



C  
C **  WRITE RESTART FILE  
C  
      IF(ISRESTO.EQ.-1.OR.ISRESTO.EQ.-11)THEN  
        CALL RESTOUT(0)  
        IF(ISTRAN(8).GE.1)THEN  
          IF(IWQRST.EQ.1) CALL WWQRST  
          IF(IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST  
        ENDIF  
      ENDIF  
      IF(ISRESTO.EQ.-2)THEN  
        CALL RESTMOD  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  COMPLETE LEAST SQUARES HARMONIC ANALYSIS  
C  
      LSLSHA=1  
      IF(ISLSHA.EQ.1) CALL LSQHARM  
C  
C**********************************************************************C  
C  
C **  OUTPUT COURANT NUMBER DIAGNOSTICS  
C  
C *** DSLLC BEGIN BLOCK
      IF(ISINWV.GT.0.AND.DEBUG)THEN
        OPEN(1,FILE='CFLMAX.OUT')  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE='CFLMAX.OUT')  
C  
        DO L=2,LA  
          WRITE(1,1991)IL(L),JL(L),(CFLUUU(L,K),K=1,KC)  
          WRITE(1,1992)(CFLVVV(L,K),K=1,KC)  
          WRITE(1,1992)(CFLWWW(L,K),K=1,KC)  
          WRITE(1,1992)(CFLCAC(L,K),K=1,KC)  
        ENDDO  
C  
        CLOSE(1)  
      ENDIF
C *** DSLLC END BLOCK
C  
 1991 FORMAT(2I5,12F8.3)  
 1992 FORMAT(10X,12F8.3)  
 1993 FORMAT(2I5,E13.5)  
C  
C**********************************************************************C  
C  
C **  OUTPUT COSMETIC VOLUME LOSSES FORM DRY CELLS  
C  
      IF(NDRYSTP.LT.0.AND.DEBUG) THEN  
C  
        OPEN(1,FILE='DRYLOSS.OUT')  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE='DRYLOSS.OUT')  
C  
        DO L=2,LA  
          WRITE(1,1993)IL(L),JL(L),VDWASTE(L)  
        ENDDO  
C  
        CLOSE(1)  
C  
      ENDIF  
C  
C**********************************************************************C  
C  
C **  OUTPUT FINAL FOOD CHAIN AVERAGING PERIOD  
C  
      IF(ISTRAN(5).GE.1.AND.ISFDCH.GE.1)CALL FOODCHAIN(1)  
C  
C**********************************************************************C  
C  
C **  OUTPUT FINAL MASS AND VOLUME BALANCES  
C  
      IF(IS2TIM.GE.1) THEN  
        IF(ISBAL.GE.1)THEN  
          CALL BAL2T5  
        ENDIF  
      ENDIF  
C  
C**********************************************************************C  
C  

      CLOSE(90)
      CLOSE(98)

      RETURN  
      END  
