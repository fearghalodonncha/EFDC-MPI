      SUBROUTINE CALEXP2T  
!  
! **  SUBROUTINE CALEXP2T CALCULATES EXPLICIT MOMENTUM EQUATION TERMS  
! **  USING A TWO TIME LEVEL SCHEME  
! CHANGE RECORD  
!  ADDED BODY FORCES FBODYFX AND FBODYFY TO EXTERNAL MOMENTUM EQUATIONS  
!  CORRECTED ORIENTATION OF MOMENTUM FLUXES FROM SINKS AND SOURCE  
!  CORRECTED 2 LAYER (KC=-2) CURVATURE ACCELERATION CORRECTION  
!  ADDED ICK2COR,CK2UUM,CK2VVM,CK2UVM,CK2UUC,CK2VVC,CK2UVC,CK2FCX,  
!  CK2FCY TO GENERALIZE TWO LAYER MOMENTUM FLUX AND CURVATURE  
!  ACCELERATION CORRECTION  
!  MODIFIED CALCULATION OF CORIOLIS-CURVATURE ACCELERATIONS AT TIDAL  
!  OPEN BOUNDARIES  
!  ADDED VIRTUAL MOMENTUM SOURCES AND SINKS FOR SUBGRID SCALE CHANNEL  
!  INTERACTIONS, INCLUDING LOCAL VARIABLES TMPVEC1,TMPVEC2,QMCSINKX,  
!  QMCSINKY,QMCSOURX,QMSOURY  
!  ADDED DRY CELL BYPASS AND CONSISTENT INITIALIZATION OF DRY VALUES  
!
!     2008-12  SANG YUK/PMC (DSLLC) CORRECTED THE EXPLICIT INTERNAL BUOYANCY FORCINGS
!  
      USE GLOBAL  
      USE OMP_LIB
	IMPLICIT NONE
!        INTEGER IINCL,JINCL,I,J
	INTEGER::L,K,LN,LS,ID,JD,KD,NWR,IU,JU,KU,LU,NS,LNW,LSE,LL
	INTEGER::LD,NMD,LHOST,LCHNU,LW,LE,LCHNV
!        DOUBLE PRECISION LST,LEND,foo
	REAL::TMPANG,WU,WV,CACSUM,CFEFF,VEAST2,VWEST2,FCORE,FCORW
	REAL::UNORT1,USOUT1,UNORT2,USOUT2,FCORN,FCORS,VTMPATU
	REAL::UTMPATV,UMAGTMP,VMAGTMP,DZICKC,DZPU,DZPV
	REAL::RCDZF,TMPVAL,WVFACT,DETH,CI11H,CI12H,CI22H,DETU
	REAL::CI11V,CI12V,CI21V,CI22V,CI21H,CI12U,CI21U,CI22U,DETV,CI11U
	REAL::UHC,UHB,VHC,VHB,UHC1,UHB1,VHC1,VHB1,UHC2,UHB2,VHC2,VHB2
	REAL::UHB1MX,UHB1MN,VHC1MX,VHC1MN,UHC1MX,UHC1MN,VHB1MX
	REAL::VHB1MN,UHB2MX,UHB2MN,VHC2MX,VHC2MN,UHC2MX,UHC2MN,VHB2MX
	REAL::VHB2MN,BOTT,QMF,QUMF,VEAST1,VWEST1
      REAL::SFXV(LCM),SFYV(LCM)
      REAL::FXVGTMP(LCM,KC),FYVGTMP(LCM,KC)

      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::DZPC
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::TMPVEC1  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::TMPVEC2  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::FUHJ  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::FVHJ  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::QMCSINKX  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::QMCSINKY  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::QMCSOURX  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::QMCSOURY  
!
      IF(.NOT.ALLOCATED(TMPVEC1))THEN
        ALLOCATE(FUHJ(LCM,KCM))  
        ALLOCATE(FVHJ(LCM,KCM))  
        ALLOCATE(QMCSINKX(LCM,KCM))  
        ALLOCATE(QMCSINKY(LCM,KCM))  
        ALLOCATE(QMCSOURX(LCM,KCM))  
        ALLOCATE(QMCSOURY(LCM,KCM))  
        ALLOCATE(TMPVEC1(KCM))  
        ALLOCATE(TMPVEC2(KCM))  
        ALLOCATE(DZPC(LCM,KCM))
        FUHJ=0.
        FVHJ=0.
        QMCSINKX=0.
        QMCSINKY=0.
        QMCSOURX=0.
        QMCSOURY=0.
        TMPVEC1=0.
        TMPVEC2=0.
        DZPC=0.
      ENDIF
!  
      IF(ISDYNSTP.EQ.0)THEN  
        DELT=DT  
      ELSE  
        DELT=DTDYN  
      ENDIF  
!  
      IF(IS2TIM.EQ.2)THEN  
        DELT=0.5*DT  
      ENDIF  
!  
      DELTI=1./DELT  
!  
      IF(N.EQ.1.AND.DEBUG)THEN  
        OPEN(1,FILE='MFLUX.DIA')  
        CLOSE(1,STATUS='DELETE')  
      ENDIF  
!  
!**********************************************************************C  
!  
! **  INITIALIZE MOMENTUM FLUXES AND CORIOLIS TERMS  
! **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS  
!  
!----------------------------------------------------------------------C  
!

!$OMP PARALLEL PRIVATE(LN,LS,LE,LW,UHC,UHB,VHC,VHB)
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
      DO L=1,LC  
        FCAXE(L)=0.0
        FCAYE(L)=0.0
        FXE(L)=0.0
        FYE(L)=0.0
        FX(L,1:KC)=0.0
        FY(L,1:KC)=0.0
      ENDDO
!$OMP END DO NOWAIT
!  
!  
!----------------------------------------------------------------------C  
!  
      IF(IS2LMC.NE.1)THEN 
        DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              LN=LNC(L)  
              LS=LSC(L)
              LE=LEAST(L)
              LW=LWEST(L)
              UHC=0.5*(UHDY(L,K)+UHDY(LS,K))  
              UHB=0.5*(UHDY(L,K)+UHDY(LE,K))  
              VHC=0.5*(VHDX(L,K)+VHDX(LW,K))  
              VHB=0.5*(VHDX(L,K)+VHDX(LN,K))  
!  
              FUHU(L,K)=MAX(UHB,0.)*U(L,  K)  ! *** CELL CENTERED 
     &                 +MIN(UHB,0.)*U(LE,K)  
              FVHU(L,K)=MAX(VHC,0.)*U(LS, K)  
     &                 +MIN(VHC,0.)*U(L,  K)
!  
              FVHV(L,K)=MAX(VHB,0.)*V(L,  K)  ! *** CELL CENTERED
     &                 +MIN(VHB,0.)*V(LN, K)  
              FUHV(L,K)=MAX(UHC,0.)*V(LW,K)  
     &                 +MIN(UHC,0.)*V(L,  K)
            ELSE
              FUHU(L,K)=0.  
              
              FVHU(L,K)=0.  
              FVHV(L,K)=0.  
              FUHV(L,K)=0.  
            ENDIF  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO

!
      END IF    !! FOR OPENMP COMPUTE BREAK INTO TWO if STATEMENTS INSTEAD OF IF - ELSEIF 
!$OMP END PARALLEL

      IF(IS2LMC.EQ.1)THEN  
!
        DO L=2,LA  
          IF(LMASKDRY(L))THEN  
            LN=LNC(L)  
            LS=LSC(L)
            LE=LEAST(L)
            LW=LWEST(L)
            UHC1=0.5*(UHDY(L,1)+UHDY(LS,1))  
            UHB1=0.5*(UHDY(L,1)+UHDY(LE,1))  
            VHC1=0.5*(VHDX(L,1)+VHDX(LW,1))  
            VHB1=0.5*(VHDX(L,1)+VHDX(LN,1))  
            UHC2=0.5*(UHDY(L,2)+UHDY(LS,2))  
            UHB2=0.5*(UHDY(L,2)+UHDY(LE,2))  
            VHC2=0.5*(VHDX(L,2)+VHDX(LW,2))  
            VHB2=0.5*(VHDX(L,2)+VHDX(LN,2))  
            UHB1MX=0.  
            UHB1MN=0.  
            VHC1MX=0.  
            VHC1MN=0.  
            UHC1MX=0.  
            UHC1MN=0.  
            VHB1MX=0.  
            VHB1MN=0.  
            UHB2MX=0.  
            UHB2MN=0.  
            VHC2MX=0.  
            VHC2MN=0.  
            UHC2MX=0.  
            UHC2MN=0.  
            VHB2MX=0.  
            VHB2MN=0.  
!  
            BOTT=ABS(UHB1*U(L,1))  
            IF(BOTT.GT.0.0)  
     &        UHB1MX=1.+CK2UUM*(UHB2-UHB1)*(U(L,2)-U(L,1))/UHB1*U(L,1)  
            BOTT=ABS(UHB1*U(LE,1))  
            IF(BOTT.GT.0.0)
     &        UHB1MN=1.+CK2UUM*(UHB2-UHB1)*(U(LE,2)-U(LE,1))/UHB1*U(LE,1)  
            BOTT=ABS(VHC1*U(LS,1))  
            IF(BOTT.GT.0.0)  
     &        VHC1MX=1.+CK2UVM*(VHC2-VHC1)*(U(LS,2)-U(LS,1))/VHC1*U(LS,1)  
            BOTT=ABS(VHC1*U(L,1))  
            IF(BOTT.GT.0.0)
     &        VHC1MN=1.+CK2UVM*(VHC2-VHC1)*(U(L,2)-U(L,1))/VHC1*U(L,1)
            BOTT=ABS(UHC1*V(LW,1))  
            IF(BOTT.GT.0.0)  
     &        UHC1MX=1.+CK2UVM*(UHC2-UHC1)*(V(LW,2)-V(LW,1))/UHC1*V(LW,1)  
            BOTT=ABS(UHC1*V(L,1))  
            IF(BOTT.GT.0.0)  
     &        UHC1MN=1.+CK2UVM*(UHC2-UHC1)*(V(L,2)-V(L,1))/UHC1*V(L,1)  
            BOTT=ABS(VHB1*V(L,1))  
            IF(BOTT.GT.0.0)  
     &        VHB1MX=1.+CK2VVM*(VHB2-VHB1)*(V(L,2)-V(L,1))/VHB1*V(L,1)  
            BOTT=ABS(VHB1*V(LN,1))  
            IF(BOTT.GT.0.0)  
     &        VHB1MN=1.+CK2VVM*(VHB2-VHB1)*(V(LN,2)-V(LN,1))/VHB1*V(LN,1)  
            BOTT=ABS(UHB2*U(L,2))
            IF(BOTT.GT.0.0)
     &        UHB2MX=1.+CK2UUM*(UHB2-UHB1)*(U(L,2)-U(L,1))/UHB2*U(L,2)  
            BOTT=ABS(UHB2*U(LE,2))  
            IF(BOTT.GT.0.0)  
     &        UHB2MN=1.+CK2UUM*(UHB2-UHB1)*(U(LE,2)-U(LE,1))/UHB2*U(LE,2)  
            BOTT=ABS(VHC2*U(LS,2))  
            IF(BOTT.GT.0.0)
     &        VHC2MX=1.+CK2UVM*(VHC2-VHC1)*(U(LS,2)-U(LS,1))/VHC2*U(LS,2)  
            BOTT=ABS(VHC2*U(L,2))  
            IF(BOTT.GT.0.0)
     &        VHC2MN=1.+CK2UVM*(VHC2-VHC1)*(U(L,2)-U(L,1))/VHC2*U(L,2)  
            BOTT=ABS(UHC2*V(LW,2))  
            IF(BOTT.GT.0.0)  
     &        UHC2MX=1.+CK2UVM*(UHC2-UHC1)*(V(LW,2)-V(LW,1))/UHC2*V(LW,2)  
            BOTT=ABS(UHC2*V(L,2))  
            IF(BOTT.GT.0.0)  
     &        UHC2MN=1.+CK2UVM*(UHC2-UHC1)*(V(L,2)-V(L,1))/UHC2*V(L,2)  
            BOTT=ABS(VHB2*V(L,2))  
            IF(BOTT.GT.0.0)  
     &        VHB2MX=1.+CK2VVM*(VHB2-VHB1)*(V(L,2)-V(L,1))/VHB2*V(L,2)  
            BOTT=ABS(VHB2*V(LN,2))  
            IF(BOTT.GT.0.0)  
     &        VHB2MN=1.+CK2VVM*(VHB2-VHB1)*(V(LN,2)-V(LN,1))/VHB2*V(LN,2)  
!  
            FUHU(L,1)=UHB1MX*MAX(UHB1,0.)*U(L ,1)  
     &               +UHB1MN*MIN(UHB1,0.)*U(LE,1)  
            FVHU(L,1)=VHC1MX*MAX(VHC1,0.)*U(LS,1)  
     &               +VHC1MN*MIN(VHC1,0.)*U(L ,1)  
            FUHV(L,1)=UHC1MX*MAX(UHC1,0.)*V(LW,1)  
     &               +UHC1MN*MIN(UHC1,0.)*V(L ,1)  
            FVHV(L,1)=VHB1MX*MAX(VHB1,0.)*V(L ,1)  
     &               +VHB1MN*MIN(VHB1,0.)*V(LN,1)  
            FUHJ(L,1)=0.  
            FVHJ(L,1)=0.  
            FUHU(L,2)=UHB2MX*MAX(UHB2,0.)*U(L ,2)  
     &               +UHB2MN*MIN(UHB2,0.)*U(LE,2)  
            FVHU(L,2)=VHC2MX*MAX(VHC2,0.)*U(LS,2)  
     &               +VHC2MN*MIN(VHC2,0.)*U(L ,2)  
            FUHV(L,2)=UHC2MX*MAX(UHC2,0.)*V(LW,2)  
     &               +UHC2MN*MIN(UHC2,0.)*V(L ,2)  
            FVHV(L,2)=VHB2MX*MAX(VHB2,0.)*V(L ,2)  
     &               +VHB2MN*MIN(VHB2,0.)*V(LN,2)  
            FUHJ(L,2)=0.  
            FVHJ(L,2)=0.  
          ENDIF  
        ENDDO  
      ENDIF  
!  
! ADD RETURN FLOW MOMENTUM FLUX  
!  
      DO NWR=1,NQWR  
        IF(NQWRMFU(NWR).GT.0)THEN  
          IU=IQWRU(NWR)  
          JU=JQWRU(NWR)  
          KU=KQWRU(NWR)  
          LU=LIJ(IU,JU)  
          NS=NQWRSERQ(NWR)  
          QMF=QWR(NWR)+QWRSERT(NS)  
          QUMF=QMF*QMF/(H1P(LU)*DZC(KU)*DZC(KU)*BQWRMFU(NWR))  
          IF(NQWRMFU(NWR).EQ.1)  FUHJ(LU     ,KU)=QUMF  
          IF(NQWRMFU(NWR).EQ.2)  FVHJ(LU     ,KU)=QUMF  
          IF(NQWRMFU(NWR).EQ.3)  FUHJ(LU+1   ,KU)=QUMF  
          IF(NQWRMFU(NWR).EQ.4)  FVHJ(LNC(LU),KU)=QUMF  
          IF(NQWRMFU(NWR).EQ.-1) FUHJ(LU     ,KU)=-QUMF  
          IF(NQWRMFU(NWR).EQ.-2) FVHJ(LU     ,KU)=-QUMF  
          IF(NQWRMFU(NWR).EQ.-3) FUHJ(LU+1   ,KU)=-QUMF  
          IF(NQWRMFU(NWR).EQ.-4) FVHJ(LNC(LU),KU)=-QUMF  
        ENDIF  
        IF(NQWRMFD(NWR).GT.0)THEN  
          ID=IQWRD(NWR)  
          JD=JQWRD(NWR)  
          KD=KQWRD(NWR)  
          LD=LIJ(ID,JD)  
          TMPANG=0.017453*ANGWRMFD(NWR)  
          TMPANG=COS(TMPANG)  
          NS=NQWRSERQ(NWR)  
          QMF=QWR(NWR)+QWRSERT(NS)  
          QUMF=TMPANG*QMF*QMF/(H1P(LD)*DZC(KD)*DZC(KD)*BQWRMFD(NWR))  
          IF(NQWRMFD(NWR).EQ.1)  FUHJ(LD     ,KD)=-QUMF  
          IF(NQWRMFD(NWR).EQ.2)  FVHJ(LD     ,KD)=-QUMF
          IF(NQWRMFD(NWR).EQ.3)  FUHJ(LD+1   ,KD)=-QUMF  
          IF(NQWRMFD(NWR).EQ.4)  FVHJ(LNC(LD),KD)=-QUMF  
          IF(NQWRMFD(NWR).EQ.-1) FUHJ(LD     ,KD)=QUMF  
          IF(NQWRMFD(NWR).EQ.-2) FVHJ(LD     ,KD)=QUMF  
          IF(NQWRMFD(NWR).EQ.-3) FUHJ(LD+1   ,KD)=QUMF  
          IF(NQWRMFD(NWR).EQ.-4) FVHJ(LNC(LD),KD)=QUMF  
!         IF(N.LE.4.AND.DEBUG)THEN  
!           WRITE(1,1112)N,NWR,NS,ID,JD,KD,NQWRMFD(NWR),H1P(LD),QMF,  
!     &                  QUMF,FUHJ(LD,KD),FVHJ(LD,KD)  
!         ENDIF  
        ENDIF  
      ENDDO  
!    
!----------------------------------------------------------------------C  
!  
! *** COMPUTE VERTICAL ACCELERATIONS
!

!$OMP PARALLEL PRIVATE(LS,WU,WV)
      DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=2,LA  
          IF(LMASKDRY(L))THEN  
            LS=LSC(L)
            LW=LWEST(L)
            WU=0.5*DXYU(L)*(W(L,K)+W(LW,K))  
            WV=0.5*DXYV(L)*(W(L,K)+W(LS,K))  

            FWU(L,K)=MAX(WU,0.)*U(L,K)  
     &          +MIN(WU,0.)*U(L,K+1)  
            FWV(L,K)=MAX(WV,0.)*V(L,K)  
     &          +MIN(WV,0.)*V(L,K+1)  
          ELSE
            FWU(L,K)=0.
            FWV(L,K)=0.
          ENDIF
  
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
!$OMP END PARALLEL

!  
!**********************************************************************C  
!  
! ** BLOCK MOMENTUM FLUX ON LAND SIDE OF TRIANGULAR CELLS  
!  
      IF(ITRICELL.GT.0)THEN
        DO K=1,KC  
          DO L=1,LA  
            FUHU(L,K)=STCUV(L)*FUHU(L,K)  
            FVHV(L,K)=STCUV(L)*FVHV(L,K)  
          ENDDO  
        ENDDO  
      ENDIF
!  
!**********************************************************************C  
!  
! **  CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS  
!  
!----------------------------------------------------------------------C  
!  
      CACSUM=0. 
      CFMAX=CF  
      IF(ISCURVATURE)THEN

        IF(ISDCCA.EQ.0)THEN  
!  
          DO K=1,KC  
            DO L=2,LA  
              IF(LMASKDRY(L))THEN  
                LN=LNC(L)  
                LE=LEAST(L)
                CAC(L,K)=( FCORC(L)*DXYP(L)  
     &            +0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L)  
     &            -0.5*SNLT*(U(LE,K)+U(L,K))*DXDJ(L) )*HP(L)  
              ELSE
                CAC(L,K)=0.0  ! *** DSLLC SINGLE LINE
              ENDIF
              CACSUM=CACSUM+CAC(L,K) 
            ENDDO  
          ENDDO  
!  
        ELSE  
!  
!  
!$OMP PARALLEL PRIVATE(LN,LE,CFEFF,CFMAX)
          DO K=1,KC  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE) REDUCTION(+:CACSUM)
            DO L=2,LA  
              LN=LNC(L)  
              LE=LEAST(L)
              CAC(L,K)=( FCORC(L)*DXYP(L)  
     &          +0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L)  
     &          -0.5*SNLT*(U(LE,K)+U(L,K))*DXDJ(L) )*HP(L)  
              CFEFF=ABS(CAC(L,K))*DXYIP(L)*HPI(L)  
              CFMAX=MAX(CFMAX,CFEFF)  
              CACSUM=CACSUM+CAC(L,K) 
            ENDDO  
!$OMP END DO NOWAIT
          ENDDO  
        ! *** ENSURE FCAY & FCAX ARE RESET
        CACSUM=ABS(CACSUM)
        IF(CACSUM.LT.1.E-7)THEN
          DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
            DO L=2,LA
              FCAX(L,K)=0.
              FCAY(L,K)=0.
            ENDDO
!$OMP END DO NOWAIT
          ENDDO
        ENDIF

!$OMP END PARALLEL

!  
          IF(N.EQ.NTS.AND.DEBUG)THEN  
            OPEN(1,FILE='CORC1.DIA')  
            CLOSE(1,STATUS='DELETE')  
            OPEN(1,FILE='CORC1.DIA')  
            K=1  
            DO L=2,LA  
              LN=LNC(L)
              LE=LEAST(L)
              WRITE(1,1111)IL(L),JL(L),LN,V(LN,K),V(L,K),DYU(LE),  
     &         DYU(L),U(LE,K),U(L,K),DXV(LN),DXV(L),HP(L),CAC(L,K)  
            ENDDO  
            CLOSE(1)  
          ENDIF  
        ENDIF  

            
      ENDIF
!  
 1111 FORMAT(3I5,10E13.4)  
 1113 FORMAT(2I5,10E13.4)  
!  
!**********************************************************************C  
!  
! **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS  
!  
!----------------------------------------------------------------------C  
!  
! **  STANDARD CALCULATION  
!   
      IF(IS2LMC.EQ.0.AND.CACSUM.GT.1.E-7)THEN
!$OMP PARALLEL PRIVATE(LN,LS,LNW,LSE,LE,LW)          
        DO K=1,KC  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE) 
          DO L=2,LA  
            IF(LMASKDRY(L))THEN
              LN=LNC(L)  
              LS=LSC(L)  
              LNW=LNWC(L)  
              LSE=LSEC(L)  
              LE=LEAST(L)
              LW=LWEST(L)
              FCAX(L,K)=0.25*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
     &            +CAC(LWEST(L),K)*(V(LNW,K)+V(LW,K)))  
              FCAY(L,K)=0.25*SCAY(L)*(CAC(L,K)*(U(LE,K)+U(L,K))  
     &            +CAC(LS,K)*(U(LSE,K)+U(LS,K)))  
            ELSE
              FCAX(L,K)=0.
              FCAY(L,K)=0.
            ENDIF  
          ENDDO  
        ENDDO  
!$OMP END PARALLEL

!  
!----------------------------------------------------------------------C  
!  
! **  MODIFICATION FOR TYPE 2 OPEN BOUNDARIES  
!  
        DO LL=1,NPBW  
          IF(ISPBW(LL).EQ.2)THEN  
            L=LPBW(LL)+1  
            LN=LNC(L)  
            DO K=1,KC  
              FCAX(L,K)=0.5*SCAX(L)*CAC(L,K)*(V(LN,K)+V(L,K))  
            ENDDO  
          ENDIF  
        ENDDO  
!  
        DO LL=1,NPBE  
          IF(ISPBE(LL).EQ.2)THEN  
            L=LPBE(LL)  
            LNW=LNWC(L)
            LW=LWEST(L)
            DO K=1,KC  
              FCAX(L,K)=0.5*SCAX(L)*CAC(LW,K)*(V(LNW,K)+V(LW,K))  
            ENDDO  
          ENDIF  
        ENDDO  
!  
        DO LL=1,NPBS  
          IF(ISPBS(LL).EQ.2)THEN  
            L=LNC(LPBS(LL))
            LE=LEAST(L)
            DO K=1,KC  
              FCAY(L,K)=0.5*SCAY(L)*CAC(L,K)*(U(LE,K)+U(L,K))  
            ENDDO  
          ENDIF  
        ENDDO  
!  
        DO LL=1,NPBN  
          IF(ISPBN(LL).EQ.2)THEN  
            L=LPBN(LL)  
            LS=LSC(L)  
            LSE=LSEC(L)  
            DO K=1,KC  
              FCAY(L,K)=0.5*SCAY(L)*CAC(LS,K)*(U(LSE,K)+U(LS,K))  
            ENDDO  
          ENDIF  
        ENDDO  
      ENDIF
!  
!----------------------------------------------------------------------C  
!  
! *** CALCULATION FOR MOMENTUM-CURVATURE CORRECTION  
! *** PMC - USED TO BE ONLY FOR 2 LAYERS, JH ALLOWED ANY # OF LAYERS
!  
      IF(IS2LMC.EQ.1.AND.CACSUM.GT.1.E-7)THEN  
!JH     IF(KC.EQ.2)THEN  
!$OMP PARALLEL DO SCHEDULE(STATIC,CHUNKSIZE) PRIVATE(LN,LS,LNW,LSE,LE,LW,
!$OMP& VEAST1,VWEST1,VEAST2,VWEST2,FCORE,FCORW,
!$OMP& UNORT1,USOUT1,UNORT2,USOUT2,FCORN,FCORS)   
        DO L=2,LA  
          IF(LMASKDRY(L))THEN  
            LN=LNC(L)  
            LS=LSC(L)  
            LNW=LNWC(L)  
            LSE=LSEC(L)  
            LW=LWEST(L)
            LE=LEAST(L)
!  
            VEAST1=V(LN,1)+V(L,1)  
            VWEST1=V(LNW,1)+V(LW,1)  
            VEAST2=V(LN,2)+V(L,2)  
            VWEST2=V(LNW,2)+V(LW,2)  
            FCORE=CK2FCX*(CAC(L,2)-CAC(L,1))*(VEAST2-VEAST1)  
            FCORW=CK2FCX*(CAC(LW,2)-CAC(LW,1))*(VWEST2-VWEST1)  
!  
            FCAX(L,1)=0.25*SCAX(L)*(  
     &                   CAC(L,1)*VEAST1+FCORE  
     &                  +CAC(LW,1)*VWEST1+FCORW)  
!  
            FCAX(L,2)=0.25*SCAX(L)*(  
     &                   CAC(L,2)*VEAST2+FCORE  
     &                  +CAC(L-2,2)*VWEST2+FCORW)  
!  
            UNORT1=U(LE,1)+U(L,1)  
            USOUT1=U(LSE,1)+U(LS,1)  
            UNORT2=U(LE,2)+U(L,2)  
            USOUT2=U(LSE,2)+U(LS,2)  
            FCORN=CK2FCY*(CAC(L,2)-CAC(L,1))*(UNORT2-UNORT1)  
            FCORS=CK2FCY*(CAC(LS,2)-CAC(LS,1))*(USOUT2-USOUT1)  
!  
            FCAY(L,1)=0.25*SCAY(L)*(  
     &                   CAC(L,1)*UNORT1+FCORN  
     &                  +CAC(LS,1)*USOUT1+FCORS)  
!  
            FCAY(L,2)=0.25*SCAY(L)*(  
     &                   CAC(L,2)*UNORT2+FCORN  
     &                  +CAC(LS,2)*USOUT2+FCORS)  
!  
          ENDIF  
        ENDDO  
      ENDIF  
!  
!----------------------------------------------------------------------C  
!  
!$OMP PARALLEL PRIVATE(LN,LS,LE,LW)   
      DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
        DO L=2,LA  
          IF(LMASKDRY(L))THEN  
            LN=LNC(L) 
            LS=LSC(L) 
            LW=LWEST(L)
            LE=LEAST(L)
            FX(L,K)=(FUHU(L,K)-FUHU(LW,K)+FVHU(LN,K)-FVHU(L,K)+FUHJ(L,K))  
            FY(L,K)=(FUHV(LE,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K)+FVHJ(L,K))
          ELSE
            FX(L,K)=0.
            FY(L,K)=0.
          ENDIF  
        ENDDO  
      ENDDO
!$OMP END PARALLEL
      ! *** TREAT BC'S NEAR EDGES
      DO LL=1,NBCS
        ! *** BC CELL
        L=LBCS(LL)
        DO K=1,KC
          FX(L,K)=SAAX(L)*FX(L,K)
          FY(L,K)=SAAY(L)*FY(L,K)
        ENDDO

        ! *** EAST/WEST ADJACENT CELL
        L=LBERC(LL)
        DO K=1,KC
          FX(L,K)=SAAX(L)*FX(L,K)
        ENDDO

        ! *** NORTH/SOUTH ADJACENT CELL
        L=LBNRC(LL)
        DO K=1,KC
          FY(L,K)=SAAY(L)*FY(L,K)
        ENDDO
      ENDDO  
!  
!----------------------------------------------------------------------C  
!  
! **  CORIOLIS-CURVATURE DIAGNOSTICS  
!  
      IF(ISDCCA.EQ.1.AND.DEBUG)THEN  
        IF(N.EQ.NTS)THEN  
          OPEN(1,FILE='CORC2.DIA')  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE='CORC2.DIA')  
          K=1  
          DO L=2,LA  
            LN=LNC(L)  
            LS=LSC(L)  
            LW=LWEST(L)
            LNW=LNWC(L)  
            LSE=LSEC(L)  
            WRITE(1,1113)IL(L),JL(L),CAC(L,K),V(LN,K),V(L,K),  
     &          CAC(LW,K),V(LNW,K),V(LW,K)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
!  
        IF(N.EQ.NTS)THEN  
          OPEN(1,FILE='CORC3.DIA')  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE='CORC3.DIA')  
          K=1  
          DO L=2,LA  
            LN=LNC(L)  
            LS=LSC(L)  
            LNW=LNWC(L)  
            LSE=LSEC(L)
            LE=LEAST(L)
            WRITE(1,1113)IL(L),JL(L),CAC(L,K),U(LE,K),U(L,K),  
     &          CAC(LS,K),U(LSE,K),U(LS,K)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
!  
        IF(N.EQ.NTS)THEN  
          OPEN(1,FILE='CORC4.DIA')  
          CLOSE(1,STATUS='DELETE')  
          OPEN(1,FILE='CORC4.DIA')  
          DO L=2,LA  
            WRITE(1,1113)IL(L),JL(L),(FCAX(L,K),K=1,KC)  
          ENDDO  
          DO L=2,LA  
            WRITE(1,1113)IL(L),JL(L),(FCAY(L,K),K=1,KC)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
      ENDIF  

!**********************************************************************C  
!  
! **  ADD VEGETATION DRAG TO HORIZONTAL ADVECTIVE ACCELERATIONS  
!  
!----------------------------------------------------------------------C   
!!!Begin SCJ block
      IF(ISVEG>=1)THEN
        FXVEGE=0.0;FYVEGE=0.0 !Zero external mode forces
        DO L=2,LA  !loop over the model area
          LW=LWEST(L)  !west cell
          LE=LEAST(L)  !east cell
          LS=LSC(L)    !south cell
          LN=LNC(L)    !north cell
          LNW=LNWC(L)  !northwest cell
          LSE=LSEC(L)  !southeast cell
          IF(.NOT.LMASKDRY(L).OR.MVEGL(L)==MVEGOW)CYCLE  !if the cell is dry, or if it is open water, or if there is no vegetation in cell L, skip this cell
          IF(MVEGL(L)==0.AND.MVEGL(LW)==0.AND.MVEGL(LS)==0)CYCLE !if not this cell and no surrounding cells are vegetation, skip
          DO K=1,KC  !loop over the model layers
            IF(RMAC(L,K)==0.0)CYCLE
            UTMPATV=0.25*(U(L,K)+U(LE,K) + U(LS,K)+U(LSE,K))  !u-velocity at v face
            VTMPATU=0.25*(V(L,K)+V(LW,K) + V(LN,K)+V(LNW,K))  !v-velocity at u face
            UMAGTMP=SQRT( U(L,K)*U(L,K)  +VTMPATU*VTMPATU )   !u-face velocity vector magnitude
            VMAGTMP=SQRT( UTMPATV*UTMPATV+  V(L,K)*V(L,K) )   !v-face velocity vector magnitude
!FXVEG/FYVEG come from CALTBXY as unitless, they are terms accounting for the drag coefficient (C_d) with a term accounting for the vegetative area density RDLSPQ (N/L^2) multiplied by the flow-facing area of each "reed" (BPVEG*HPVEG*RDLPSQ)
!FXVEG/FYVEG only change inasmuch as the water depth changes and are zero in layers not penetrated by vegetation. Layer thickness is included by multiplying by DZC(K) below. 
!FXVEG/FYVEG are C_d(N/L^2)*BPVEG*HPVEG (plus some other factors) [-], this dimensionless term is now called C_v
!FXVEG/FYVEG are now multiplied by the face-centered cell area and local cell-averaged velocity to yield C_v|q|A in units of [m^3/s]
            FXVEG (L,K)=  UMAGTMP * SUB(L) * FXVEG(L,K)*DZC(K) * DXYU(L) ![m^3/s] |q_x|C_d(N/L^2)A Note that FXVEG is multiplied by DZC to divide out by the HP term in CALTBXY (actual layer force)
            FYVEG (L,K)=  VMAGTMP * SVB(L) * FYVEG(L,K)*DZC(K) * DXYV(L) ![m^3/s] |q_y|C_d(N/L^2)A Note that FYVEG is multiplied by DZC to divide out by the HP term in CALTBXY (actual layer force)
            FXVEGE(L)  =FXVEGE(L)+FXVEG(L,K)*HUI(L)*DXYIU(L) ![1/s] Calculate vegetative dissipation for FUHDYE for momentum conservation in CALPUV (need to have sum of forces in column, not average, provided to CALPUV). This is volume normalized as CALPUV is expecting.
            FYVEGE(L)  =FYVEGE(L)+FYVEG(L,K)*HVI(L)*DXYIV(L) ![1/s] Calculate vegetative dissipation for FVHDXE for momentum conservation in CALPUV (need to have sum of forces in column, not average, provided to CALPUV). This is volume normalized as CALPUV is expecting.
!Note that FXVEG/FYVEG come from CALTBXY as the drag coefficient that has been multiplied by the water-column height
          ENDDO
!FXVEG/FYVEG are added to the body forces as C_d(N/L^2)A|q|q
!FXVEG/FYVEG are multiplied by the local velocity to yield units of [m^4/s^2]
          SFXV(L)=SUM(FXVEG(L,1:KC)*DZC(1:KC))!*FLOAT(KC)
          SFYV(L)=SUM(FYVEG(L,1:KC)*DZC(1:KC))!*FLOAT(KC)
          FXVGTMP(L,1:KC)=(FXVEG(L,1:KC)-SFXV(L))*U(L,1:KC)
          FYVGTMP(L,1:KC)=(FYVEG(L,1:KC)-SFYV(L))*V(L,1:KC)
          FX(L,1:KC)=FX(L,1:KC)+FXVGTMP(L,1:KC) ![m^4/s^2] adding vegetative resistance to the body force (no net force added) to induce layer shear. 
          FY(L,1:KC)=FY(L,1:KC)+FYVGTMP(L,1:KC) ![m^4/s^2] adding vegetative resistance to the body force (no net force added) to induce layer shear.
        ENDDO
        IF(LMHK)THEN
          CALL MHKPWRDIS !MHK devices exist
          FXVEGE(:)=FXVEGE(:)+FXMHKE(:)+FXSUPE(:) !Add MHK to vegetative dissipation in FUHDYE for momentum conservation in CALPUV (MHK/SUP come in from MHKPWR as volume normalized)
          FYVEGE(:)=FYVEGE(:)+FYMHKE(:)+FYSUPE(:) !Add MHK to vegetative dissipation in FVHDXE for momentum conservation in CALPUV (MHK/SUP come in from MHKPWR as volume normalized)
        ENDIF
      ENDIF
!!!End SCJ block
 1947 FORMAT(3I5,10E12.4)
 1948 FORMAT(15X,10E12.4)

! FEARGHAL aquaculture installation simulation capabilities
      IF (aquaincl == 1) CALL AQUAEXTR

!**********************************************************************C  
!  
! **  ADD HORIZONTAL MOMENTUM DIFFUSION TO ADVECTIVE ACCELERATIONS  
!  
!**********************************************************************C  
       IF(ISHDMF.GE.1)THEN 
!$OMP PARALLEL
        DO K=1,KC  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            IF(LMASKDRY(L))THEN
              FX(L,K)=FX(L,K)-(FMDUX(L,K)+FMDUY(L,K))  
              FY(L,K)=FY(L,K)-(FMDVX(L,K)+FMDVY(L,K))  
            ENDIF
          ENDDO  
        ENDDO
!$OMP END PARALLEL
      ENDIF  

!**********************************************************************C  
!  
! **  ADD BODY FORCE TO ADVECTIVE ACCELERATIONS  
! **  DISTRIBUTE UNIFORMLY OVER ALL LAYERS IF ISBODYF=1  
! **  DISTRIBUTE OVER SURFACE LAYER IF ISBODYF=2  
!  
!----------------------------------------------------------------------C  
      IF(ISBODYF.EQ.1.OR.ISUVDA.GE.1)THEN  
!  
        DO K=1,KC  
          DO L=2,LA  
            FX(L,K)=FX(L,K)-DYU(L)*HU(L)*FBODYFX(L,K)  
            FY(L,K)=FY(L,K)-DXV(L)*HV(L)*FBODYFY(L,K)  
          ENDDO  
        ENDDO  
!  
      ENDIF  
!  
      IF(ISBODYF.EQ.2)THEN  
!  
        DZICKC=1./DZC(KC)  
        DO L=2,LA  
          FX(L,KC)=FX(L,KC)-DZICKC*DYU(L)*HU(L)*FBODYFX(L,K)  
          FY(L,KC)=FY(L,KC)-DZICKC*DXV(L)*HV(L)*FBODYFY(L,K)  
        ENDDO  
!  
      ENDIF  
!  
!**********************************************************************C  
!  
! ** ADD EXPLICIT NONHYDROSTATIC PRESSURE  
!  
      IF(KC.GT.1.AND.ISPNHYDS.GE.1) THEN  
!  
        TMPVAL=2./(DZC(1)+DZC(2))  
        DO L=2,LA  
          DZPC(L,1)=TMPVAL*(PNHYDS(L,2)-PNHYDS(L,1))  
        ENDDO  
!  
        TMPVAL=2./(DZC(KC)+DZC(KC-1))  
        DO L=2,LA  
          DZPC(L,KC)=TMPVAL*(PNHYDS(L,KC)-PNHYDS(L,KC-1))  
        ENDDO  
  
        IF(KC.GE.3)THEN  
          DO K=2,KS  
            TMPVAL=2./(DZC(K+1)+2.*DZC(K)+DZC(K-1))  
            DO L=2,LA  
              DZPC(L,K)=TMPVAL*(PNHYDS(L,K+1)-PNHYDS(L,K-1))  
            ENDDO  
          ENDDO  
        ENDIF  
!  
        DO K=1,KC  
          DO L=2,LA  
            LS=LSC(L)  
            DZPU=0.5*(DZPC(L,K)+DZPC(LWEST(L),K))  
            DZPV=0.5*(DZPC(L,K)+DZPC(LS ,K))  
            FX(L,K)=FX(L,K)+SUB(L)*DYU(L)*  
     &          ( HU(L)*(PNHYDS(L,K)-PNHYDS(LWEST(L),K))  
     &          -( BELV(L)-BELV(LWEST(L))+ZZ(K)*(HP(L)-HP(LWEST(L))) )*DZPU )  
            FY(L,K)=FY(L,K)+SVB(L)*DXV(L)*  
     &          ( HV(L)*(PNHYDS(L,K)-PNHYDS(LS ,K))  
     &          -( BELV(L)-BELV(LS )+ZZ(K)*(HP(L)-HP(LS )) )*DZPV )  
          ENDDO  
        ENDDO  
!  
      ENDIF  
C  
C----------------------------------------------------------------------C  
C  
C **  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.  
C  
C *** DSLLC BEGIN BLOCK
      IF(ISWAVE.EQ.2)THEN
C
        IF(N.LT.NTSWV)THEN  
          TMPVAL=FLOAT(N)/FLOAT(NTSWV)  
          WVFACT=0.5-0.5*COS(PI*TMPVAL)  
        ELSE  
          WVFACT=1.0  
        ENDIF  
C
        IF(ISDRY.GT.0)THEN
          DO K=1,KC
            DO L=2,LA
              IF(LMASKDRY(L))THEN  
                FX(L,K)=FX(L,K)+WVFACT*SAAX(L)*FXWAVE(L,K)
                FY(L,K)=FY(L,K)+WVFACT*SAAY(L)*FYWAVE(L,K)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO K=1,KC
            DO L=2,LA
              FX(L,K)=FX(L,K)+WVFACT*SAAX(L)*FXWAVE(L,K)
              FY(L,K)=FY(L,K)+WVFACT*SAAY(L)*FYWAVE(L,K)
            ENDDO
          ENDDO
        ENDIF
C  
      ENDIF  
C *** DSLLC END BLOCK
C  
C**********************************************************************C  
C  
C **  CALCULATE EXTERNAL ACCELERATIONS  
C  
C----------------------------------------------------------------------C  
C  
!$OMP PARALLEL
      IF(ISDRY.GT.0)THEN
        DO K=1,KC 
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE) 
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              FCAXE(L)=FCAXE(L)+FCAX(L,K)*DZC(K)  
              FCAYE(L)=FCAYE(L)+FCAY(L,K)*DZC(K)  
              FXE(L)=FXE(L)+FX(L,K)*DZC(K)  
              FYE(L)=FYE(L)+FY(L,K)*DZC(K)  
            ENDIF  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
      ELSE
        DO K=1,KC  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            FCAXE(L)=FCAXE(L)+FCAX(L,K)*DZC(K)  
            FCAYE(L)=FCAYE(L)+FCAY(L,K)*DZC(K)  
            FXE(L)=FXE(L)+FX(L,K)*DZC(K)  
            FYE(L)=FYE(L)+FY(L,K)*DZC(K)  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO
      ENDIF
!  
!**********************************************************************C  
!  
! **  COMPLETE CALCULATION OF INTERNAL ADVECTIVE ACCELERATIONS  
!  
!----------------------------------------------------------------------C  
!  
      IF(KC.GT.1)THEN
        DO K=1,KC  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              FX(L,K)=FX(L,K)+SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
              FY(L,K)=FY(L,K)+SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)  
            ENDIF  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
      ENDIF
!$OMP END PARALLEL

!  
!**********************************************************************C  
!  
! **  ADD SUBGRID SCALE CHANNEL VIRTURAL MOMENTUM SOURCES AND SINKS  
!  
!----------------------------------------------------------------------C  
!  
      IF(MDCHH.GE.1.AND.ISCHAN.EQ.3)THEN  
!  
        DO K=1,KC  
          DO L=2,LA  
            QMCSOURX(L,K)=0.  
            QMCSOURY(L,K)=0.  
            QMCSINKX(L,K)=0.  
            QMCSINKY(L,K)=0.  
          ENDDO  
        ENDDO  
!  
        DO NMD=1,MDCHH  
!  
          LHOST=LMDCHH(NMD)  
          LCHNU=LMDCHU(NMD)  
          LCHNV=LMDCHV(NMD)  
!  
          DETH=CUE(LHOST)*CVN(LHOST)-CUN(LHOST)*CVE(LHOST)  
          CI11H=CVN(LHOST)/DETH  
          CI12H=-CUN(LHOST)/DETH  
          CI21H=-CVE(LHOST)/DETH  
          CI22H=CUE(LHOST)/DETH  
!  
          DETU=CUE(LCHNU)*CVN(LCHNU)-CUN(LCHNU)*CVE(LCHNU)  
          CI11U=CVN(LCHNU)/DETU  
          CI12U=-CUN(LCHNU)/DETU  
          CI21U=-CVE(LCHNU)/DETU  
          CI22U=CUE(LCHNU)/DETU  
!  
          DETV=CUE(LCHNV)*CVN(LCHNV)-CUN(LCHNV)*CVE(LCHNV)  
          CI11V=CVN(LCHNV)/DETV  
          CI12V=-CUN(LCHNV)/DETV  
          CI21V=-CVE(LCHNV)/DETV  
          CI22V=CUE(LCHNV)/DETV  
!  
!         X-DIRECTION CHANNEL  
          IF(MDCHTYP(NMD).EQ.1)THEN  
            IF(QCHANU(NMD).GT.0.0)THEN  
              DO K=1,KC  
                QMCSINKX(LCHNU,K)=QMCSINKX(LCHNU,K)  
     &              -0.5*DZC(K)*QCHANU(NMD)*(U(LCHNU,K)+U(LCHNU+1,K))  
                QMCSINKY(LCHNU,K)=QMCSINKY(LCHNU,K)  
     &              -0.5*DZC(K)*QCHANU(NMD)*(V(LCHNU,K)+V(LNC(LCHNU),K))  
              ENDDO  
              DO K=1,KC  
                TMPVEC1(K)=CUE(LCHNU)*QMCSINKX(LCHNU,K)  
     &              +CVE(LCHNU)*QMCSINKY(LCHNU,K)  
                TMPVEC2(K)=CUN(LCHNU)*QMCSINKX(LCHNU,K)  
     &              +CVN(LCHNU)*QMCSINKY(LCHNU,K)  
              ENDDO  
              DO K=1,KC  
                QMCSOURX(LHOST,K)=QMCSOURX(LHOST,K)  
     &              +CI11H*TMPVEC1(K)+CI12H*TMPVEC2(K)  
                QMCSOURY(LHOST,K)=QMCSOURY(LHOST,K)  
     &              +CI21H*TMPVEC1(K)+CI22H*TMPVEC2(K)  
              ENDDO  
            ELSE  
              DO K=1,KC  
                QMCSINKX(LHOST,K)=QMCSINKX(LHOST,K)  
     &              +0.5*DZC(K)*QCHANU(NMD)*(U(LHOST,K)+U(LHOST+1,K))  
                QMCSINKY(LHOST,K)=QMCSINKY(LCHNU,K)  
     &              +0.5*DZC(K)*QCHANU(NMD)*(V(LHOST,K)+V(LNC(LHOST),K))  
              ENDDO  
              DO K=1,KC  
                TMPVEC1(K)=CUE(LHOST)*QMCSINKX(LHOST,K)  
     &              +CVE(LHOST)*QMCSINKY(LHOST,K)  
                TMPVEC2(K)=CUN(LHOST)*QMCSINKX(LCHNU,K)  
     &              +CVN(LHOST)*QMCSINKY(LHOST,K)  
              ENDDO  
              DO K=1,KC  
                QMCSOURX(LCHNU,K)=QMCSOURX(LCHNU,K)  
     &              -CI11U*TMPVEC1(K)-CI12U*TMPVEC2(K)  
                QMCSOURY(LCHNU,K)=QMCSOURY(LCHNU,K)  
     &              -CI21U*TMPVEC1(K)-CI22U*TMPVEC2(K)  
              ENDDO  
            ENDIF  
          ENDIF  
!  
!         Y-DIRECTION CHANNEL  
          IF(MDCHTYP(NMD).EQ.2)THEN  
            IF(QCHANV(NMD).GT.0.0)THEN  
              DO K=1,KC  
                QMCSINKX(LCHNV,K)=QMCSINKX(LCHNV,K)  
     &              -0.5*DZC(K)*QCHANV(NMD)*(U(LCHNV,K)+U(LCHNV+1,K))  
                QMCSINKY(LCHNV,K)=QMCSINKY(LCHNV,K)  
     &              -0.5*DZC(K)*QCHANV(NMD)*(V(LCHNV,K)+V(LNC(LCHNV),K))  
              ENDDO  
              DO K=1,KC  
                TMPVEC1(K)=CUE(LCHNV)*QMCSINKX(LCHNV,K)  
     &              +CVE(LCHNV)*QMCSINKY(LCHNV,K)  
                TMPVEC2(K)=CUN(LCHNV)*QMCSINKX(LCHNV,K)  
     &              +CVN(LCHNV)*QMCSINKY(LCHNV,K)  
              ENDDO  
              DO K=1,KC  
                QMCSOURX(LHOST,K)=QMCSOURX(LHOST,K)  
     &              +CI11H*TMPVEC1(K)+CI12H*TMPVEC2(K)  
                QMCSOURY(LHOST,K)=QMCSOURY(LHOST,K)  
     &              +CI21H*TMPVEC1(K)+CI22H*TMPVEC2(K)  
              ENDDO  
            ELSE  
              DO K=1,KC  
                QMCSINKX(LHOST,K)=QMCSINKX(LHOST,K)  
     &              +0.5*DZC(K)*QCHANV(NMD)*(U(LHOST,K)+U(LHOST+1,K))  
                QMCSINKY(LHOST,K)=QMCSINKY(LCHNV,K)  
     &              +0.5*DZC(K)*QCHANV(NMD)*(V(LHOST,K)+V(LNC(LHOST),K))  
              ENDDO  
              DO K=1,KC  
                TMPVEC1(K)=CUE(LHOST)*QMCSINKX(LHOST,K)  
     &              +CVE(LHOST)*QMCSINKY(LHOST,K)  
                TMPVEC2(K)=CUN(LHOST)*QMCSINKX(LCHNU,K)  
     &              +CVN(LHOST)*QMCSINKY(LHOST,K)  
              ENDDO  
              DO K=1,KC  
                QMCSOURX(LCHNV,K)=QMCSOURX(LCHNV,K)  
     &              -CI11V*TMPVEC1(K)-CI12V*TMPVEC2(K)  
                QMCSOURY(LCHNV,K)=QMCSOURY(LCHNV,K)  
     &              -CI21V*TMPVEC1(K)-CI22V*TMPVEC2(K)  
              ENDDO  
            ENDIF  
          ENDIF  
!  
        ENDDO  
!  
        DO K=1,KC  
          DO L=2,LA  
            IF(QMCSOURX(L,K).NE.0.0)THEN 
              LE=LEAST(L) 
              TMPVAL=SUB(L)+SUB(LE)  
              TMPVAL=MAX(TMPVAL,1.0)  
              FX(L,K)=FX(L,K)-SUB(L)*QMCSOURX(L,K)/TMPVAL  
              FX(LE,K)=FX(LE,K)-SUB(LE)*QMCSOURX(L,K)/TMPVAL  
            ENDIF  
            IF(QMCSOURY(L,K).NE.0.0)THEN  
              LN=LNC(L)  
              TMPVAL=SVB(L)+SVB(LN)  
              TMPVAL=MAX(TMPVAL,1.0)  
              FY(L,K)=FY(L,K)-SVB(L)*QMCSOURX(L,K)/TMPVAL  
              FY(LN,K)=FY(LN,K)-SVB(LN)*QMCSOURX(L,K)/TMPVAL  
            ENDIF  
            IF(QMCSINKX(L,K).NE.0.0)THEN 
              LE=LEAST(L) 
              TMPVAL=SUB(L)+SUB(LE)  
              TMPVAL=MAX(TMPVAL,1.0)  
              FX(L,K)=FX(L,K)-SUB(L)*QMCSINKX(L,K)/TMPVAL  
              FX(LE,K)=FX(LE,K)-SUB(LE)*QMCSINKX(L,K)/TMPVAL  
            ENDIF  
            IF(QMCSINKY(L,K).NE.0.0)THEN  
              LN=LNC(L)  
              TMPVAL=SVB(L)+SVB(LNC(L))  
              TMPVAL=MAX(TMPVAL,1.0)  
              FY(L,K)=FY(L,K)-SVB(L)*QMCSINKX(L,K)/TMPVAL  
              FY(LN,K)=FY(LN,K)-SVB(LN)*QMCSINKX(L,K)/TMPVAL  
            ENDIF  
          ENDDO  
        ENDDO  
!  
      ENDIF  
!  
!**********************************************************************C  
!  
! **  CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR  
! **  THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP  
! **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV  
!  
!----------------------------------------------------------------------C  
!  
c      IINTPG=0  
!  
!     ORIGINAL  
!  
      IF(BSC.GT.1.E-6)THEN
!$OMP PARALLEL PRIVATE(LS,LW)     
        IF(IINTPG.EQ.0)THEN  
          DO K=1,KS
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
            DO L=2,LA  
              LS=LSC(L)  
              LW=LWEST(L)
              FBBX(L,K)=SBX(L)*GP*HU(L)*  
     &          ( HU(L)*( (B(L,K+1)-B(LW,K+1))*DZC(K+1)  
     &          +(B(L,K)-B(LW,K))*DZC(K) )  
     &          -(B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*  
     &          (BELV(L)-BELV(LW)+Z(K)*(HP(L)-HP(LW))) )  
              FBBY(L,K)=SBY(L)*GP*HV(L)*  
     &          ( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)  
     &          +(B(L,K)-B(LS,K))*DZC(K) )  
     &          -(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*  
     &          (BELV(L)-BELV(LS)+Z(K)*(HP(L)-HP(LS))) )  
            ENDDO  
!$OMP END DO NOWAIT
          ENDDO  
!  
        ENDIF  
!  
! *** JACOBIAN
!  
        IF(IINTPG.EQ.1.)THEN
        K=1  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
        DO L=2,LA  
          LS=LSC(L)  
          LW=LWEST(L)
          FBBX(L,K)=SBX(L)*GP*HU(L)*  
     &        ( 0.5*HU(L)*( (B(L,K+2)-B(LW,K+2))*DZC(K+2)  
     &        +(B(L,K+1)-B(LW,K+1))*DZC(K+1)  
     &        +(B(L,K  )-B(LW,K  ))*DZC(K  )  
     &        +(B(L,K  )-B(LW,K  ))*DZC(K  ) )  
     &        -0.5*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))*  
     &        (BELV(L)-BELV(LW)+Z(K+1)*(HP(L)-HP(LW)))  
     &        -0.5*(B(L,K  )-B(L,K  )+B(LW,K  )-B(LW,K  ))*  
     &        (BELV(L)-BELV(LW)+Z(K-1)*(HP(L)-HP(LW))) )
!  
          FBBY(L,K)=SBY(L)*GP*HV(L)*  
     &        ( 0.5*HV(L)*( (B(L,K+2)-B(LS ,K+2))*DZC(K+2)  
     &        +(B(L,K+1)-B(LS ,K+1))*DZC(K+1)  
     &        +(B(L,K  )-B(LS ,K  ))*DZC(K  )  
     &        +(B(L,K  )-B(LS ,K  ))*DZC(K  ) )  
     &        -0.5*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*  
     &        (BELV(L)-BELV(LS)+Z(K+1)*(HP(L)-HP(LS)))  
     &        -0.5*(B(L,K  )-B(L,K  )+B(LS ,K  )-B(LS ,K  ))*  
     &        (BELV(L)-BELV(LS )+Z(K-1)*(HP(L)-HP(LS ))) )  
        ENDDO  
!$OMP END DO NOWAIT
!  
        IF(KC.GT.2)THEN  
          K=KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
          DO L=2,LA  
            LS=LSC(L)  
            LW=LWEST(L)
            FBBX(L,K)=SBX(L)*GP*HU(L)*  
     &          ( 0.5*HU(L)*( (B(L,K+1)-B(LW,K+1))*DZC(K+1)  
     &          +(B(L,K+1)-B(LW,K+1))*DZC(K+1)  
     &          +(B(L,K  )-B(LW,K  ))*DZC(K  )  
     &          +(B(L,K-1)-B(LW,K-1))*DZC(K-1) )  
     &          -0.5*(B(L,K+1)-B(L,K+1)+B(LW,K+1)-B(LW,K+1))*  
     &          (BELV(L)-BELV(LW)+Z(K+1)*(HP(L)-HP(LW)))  
     &          -0.5*(B(L,K  )-B(L,K-1)+B(LW,K  )-B(LW,K-1))*  
     &          (BELV(L)-BELV(LW)+Z(K-1)*(HP(L)-HP(LW))) )  
            FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*  
     &          ( 0.5*HV(L)*( (B(L,K+1)-B(LS ,K+1))*DZC(K+1)  
     &          +(B(L,K+1)-B(LS ,K+1))*DZC(K+1)  
     &          +(B(L,K  )-B(LS ,K  ))*DZC(K  )  
     &          +(B(L,K-1)-B(LS ,K-1))*DZC(K-1) )  
     &          -0.5*(B(L,K+1)-B(L,K+1)+B(LS ,K+1)-B(LS ,K+1))*  
     &          (BELV(L)-BELV(LS)+Z(K+1)*(HP(L)-HP(LS)))  
     &          -0.5*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))*  
     &          (BELV(L)-BELV(LS )+Z(K-1)*(HP(L)-HP(LS ))) )  
          ENDDO  
!$OMP END DO NOWAIT
        ENDIF  
!  
        IF(KC.GT.3)THEN  
          DO K=1,KS
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)    
            DO L=2,LA  
              LS=LSC(L)  
              LW=LWEST(L)
              FBBX(L,K)=SBX(L)*GP*HU(L)*  
     &            ( 0.5*HU(L)*( (B(L,K+2)-B(LW,K+2))*DZC(K+2)  
     &            +(B(L,K+1)-B(LW,K+1))*DZC(K+1)  
     &            +(B(L,K  )-B(LW,K  ))*DZC(K  )  
     &            +(B(L,K-1)-B(LW,K-1))*DZC(K-1) )  
     &            -0.5*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))*  
     &            (BELV(L)-BELV(LW)+Z(K+1)*(HP(L)-HP(LW)))  
     &            -0.5*(B(L,K  )-B(L,K-1)+B(LW,K  )-B(LW,K-1))*  
     &            (BELV(L)-BELV(LW)+Z(K-1)*(HP(L)-HP(LW))) )  
              FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*  
     &            ( 0.5*HV(L)*( (B(L,K+2)-B(LS ,K+2))*DZC(K+2)  
     &            +(B(L,K+1)-B(LS ,K+1))*DZC(K+1)  
     &            +(B(L,K  )-B(LS ,K  ))*DZC(K  )  
     &            +(B(L,K-1)-B(LS ,K-1))*DZC(K-1) )  
     &            -0.5*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*  
     &            (BELV(L)-BELV(LS)+Z(K+1)*(HP(L)-HP(LS)))  
     &            -0.5*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))*  
     &            (BELV(L)-BELV(LS )+Z(K-1)*(HP(L)-HP(LS ))) )  
            ENDDO  
!$OMP END DO NOWAIT
          ENDDO  
        ENDIF  
!  
      ENDIF  
!  
!     FINITE VOLUME  
!  
        IF(IINTPG.EQ.2)THEN  
!  
        DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)  
          DO L=2,LA  
            LS=LSC(L)  
            LW=LWEST(L)
            FBBX(L,K)=SBX(L)*GP*HU(L)*  
     &          ( ( HP(L)*B(L,K+1)-HP(LW)*B(LW,K+1) )*DZC(K+1)  
     &          +( HP(L)*B(L,K  )-HP(LW)*B(LW,K  ) )*DZC(K  ) )  
     &          -RNEW*SBX(L)*GP*(BELV(L)-BELV(LW))*  
     &          ( HP(L)*B(L,K+1)-HP(L)*B(L,K)  
     &          +HP(LW)*B(LW,K+1)-HP(LW)*B(LW,K) )  
     &          -RNEW*SBX(L)*GP*(HP(L)-HP(LW))*  
     &          ( HP(L)*ZZ(K+1)*B(L,K+1)-HP(L)*ZZ(K)*B(L,K)  
     &          +HP(LW)*ZZ(K+1)*B(LW,K+1)-HP(LW)*ZZ(K)*B(LW,K) )  
            FBBY(L,K)=SBY(L)*GP*HV(L)*  
     &          ( ( HP(L)*B(L,K+1)-HP(LS )*B(LS ,K+1) )*DZC(K+1)  
     &          +( HP(L)*B(L,K  )-HP(LS )*B(LS ,K  ) )*DZC(K  ) )  
     &          -RNEW*SBY(L)*GP*(BELV(L)-BELV(LS ))*  
     &          ( HP(L)*B(L,K+1)-HP(L)*B(L,K)  
     &          +HP(LS)*B(LS ,K+1)-HP(LS)*B(LS ,K) )  
     &          -RNEW*SBY(L)*GP*(HP(L)-HP(LS ))*  
     &          ( HP(L)*ZZ(K+1)*B(L,K+1)-HP(L)*ZZ(K)*B(L,K)  
     &          +HP(LS)*ZZ(K+1)*B(LS ,K+1)-HP(LS)*ZZ(K)*B(LS ,K) )  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
!  
        ENDIF  
!$OMP END PARALLEL 

      ENDIF  ! *** END OF BOUYANCY
!
!     IF(N.EQ.1)THEN
!       OPEN(1,FILE='BUOY.DIA',STATUS='UNKNOWN')
!       DO L=2,LA
!        DO K=1,KS
!        TMP3D(K)=SUBO(L)*FBBX(L,K)
!        ENDDO
!       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
!        DO K=1,KS
!        TMP3D(K)=SVBO(L)*FBBY(L,K)
!        ENDDO
!       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
!       ENDDO
!       CLOSE(1)
!     ENDIF
!
! 1111 FORMAT(2I5,2X,8E12.4)       
!
!**********************************************************************C
!
! **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
!
!----------------------------------------------------------------------C
!
!$OMP PARALLEL PRIVATE(RCDZF)
      IF(KC.GT.1)THEN
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=1,LC
          DU(L,KC)=0.0
          DV(L,KC)=0.0
        ENDDO  
!$OMP END DO NOWAIT
        DO K=1,KS
          RCDZF=CDZF(K)
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA
            IF(LMASKDRY(L))THEN
              DU(L,K)=RCDZF*( HU(L)*(U(L,K+1)-U(L,K))*DELTI
     &           +DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)+FBBX(L,K)
     &           +SNLT*(FX(L,K)-FX(L,K+1))) )
              DV(L,K)=RCDZF*( HV(L)*(V(L,K+1)-V(L,K))*DELTI
     &           +DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K)
     &           +SNLT*(FY(L,K)-FY(L,K+1))) )
            ELSE
              ! *** TEMPORARY VARIABLE, SO MUST BE INITIALIZED
              DU(L,K)=0.0  
              DV(L,K)=0.0  
            ENDIF
          ENDDO
!$OMP END DO NOWAIT
        ENDDO
      ENDIF
!
!      IF(ISTL.EQ.2)THEN
! 
      IF(NWSER.GT.0)THEN
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=2,LA
          DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
          DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
        ENDDO
!$OMP END DO NOWAIT
      ENDIF
!$OMP END PARALLEL

!
!      ENDIF
!
!**********************************************************************C
!
!      IF(N.LE.4)THEN
!        CLOSE(1)
!      ENDIF
!
 1112 FORMAT('N,NW,NS,I,J,K,NF,H,Q,QU,FUU,FVV=',/,2X,7I5,5E12.4)
!
!**********************************************************************C
!
      RETURN
      END