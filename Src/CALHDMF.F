      SUBROUTINE CALHDMF     
C
C *** CALDMF CALCULATES THE HORIZONTAL VISCOSITY AND     
C *** DIFFUSIVE MOMENTUM FLUXES. THE VISCOSITY, AH IS CALCULATED USING 
C *** SMAGORINSKY'S SUBGRID SCALE FORMULATION PLUS A CONSTANT AHO      
C
C *** ONLY VALID FOR ISHDMF.GE.1    
C
C CHANGE RECORD
C     REWRITTEN BY PAUL M. CRAIG  NOV/DEC 2004
C     2008-10  SANG YUK (DSLLC) CORRECTED THE DIFFUSIVE MOMENTUM FLUXES COMPUTATION
C
      USE GLOBAL      
!$      USE OMP_LIB
      USE PARALLEL_MPI
      IMPLICIT NONE
      INTEGER::L,LN,LS,LW,LE,K,LL,J,I    
      REAL::SLIPCO,DY2DZBR,DX2DZBR,CSDRAG,SLIPFAC,TMPVAL,DSQR,WVFACT   
      REAL::DTMPH,DTMPX,AHWVX,SXYLN,SXYEE     
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::AHEE  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::AHNN  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::SXY   
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::SXY2CC
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::SXY2EE
      REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)::SXY2NN
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::HMC     
      IF(.NOT.ALLOCATED(AHEE))THEN  
       ALLOCATE(AHEE(LCM,KCM))     
       ALLOCATE(AHNN(LCM,KCM))     
       ALLOCATE(SXY(LCM,KCM))      
       ALLOCATE(SXY2CC(LCM,KCM))   
       ALLOCATE(SXY2EE(LCM,KCM))   
       ALLOCATE(SXY2NN(LCM,KCM))   
       ALLOCATE(HMC(LCM))   
       AHEE=0.0      
       AHNN=0.0      
       SXY=0.0
       SXY2CC=0.0    
       SXY2EE=0.0    
       SXY2NN=0.0    
       HMC=0.0
      ENDIF    
C
      AHMAX=AHO
C
C **  CALCUATE TYPE FLAGS    
C
      IF(ISDRY.GE.1.OR.N.LT.5)THEN      
        ! *** ICORDYU  
        DO L=2,LA  
          LS=LSC(L)  
          IF(SUB(L).LT.0.5.AND.SUB(LS).LT.0.5) ICORDYU(L)=0  
          IF(SUB(L).GT.0.5.AND.SUB(LS).GT.0.5) ICORDYU(L)=1  
          IF(SUB(L).LT.0.5.AND.SUB(LS).GT.0.5) ICORDYU(L)=2  
          IF(SUB(L).GT.0.5.AND.SUB(LS).LT.0.5) ICORDYU(L)=3  
        ENDDO  

        ! *** ICORDXV  
        DO L=2,LA  
          LW=LWEST(L)  
          IF(SVB(L).LT.0.5.AND.SVB(LW).LT.0.5) ICORDXV(L)=0  
          IF(SVB(L).GT.0.5.AND.SVB(LW).GT.0.5)THEN  
            ICORDXV(L)=1  
            IF(SUB(L).LT.0.5) ICORDXV(L)=3  
          ENDIF  
          IF(SVB(L).LT.0.5.AND.SVB(LW).GT.0.5) ICORDXV(L)=2  
          IF(SVB(L).GT.0.5.AND.SVB(LW).LT.0.5) ICORDXV(L)=3  
        ENDDO  
      ENDIF
C  
C **  CALCULATE HORIZONTAL VELOCITY SHEARS
C
#ifdef key_mpi
      CALL COMMUNICATE_3D(U)
      CALL COMMUNICATE_3D(V)
#endif
      ! *** SXX+SYY DEFINED AT CELL CENTERS AND STORED IN DXU1(L,K)    
      IF(AHD.GT.0.0)THEN     
        SLIPCO=0.5/SQRT(AHD) 
      ENDIF    
C
      ! *** DYU1 = dU/dY
!$OMP PARALLEL PRIVATE(LS,LN,LW,LE,DY2DZBR,CSDRAG,SLIPFAC)
      DO K=1,KC  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=2,LA  
          LS=LSC(L)  
          LN=LNC(L)   
          LW=LWEST(L)  
          LE=LEAST(L)
          ! *** DXU1 = dU/dX, UNITS: 1/S   
          DXU1(L,K)=SUB(L)*SUB(LE )*(U(LE ,K)-U(L,K))/DXP(L)    
          ! *** DYV1 = dV/dY, UNITS: 1/S   
          DYV1(L,K)=SVB(L)*SVB(LN )*(V(LN ,K)-V(L,K))/DYP(L)     
          IF(ICORDYU(L).EQ.1)THEN
           DYU1(L,K)=2.*SVB(L)*SVB(LS)*(U(L,K)-U(LS,K))/(DYU(L)+DYU(LS))
          ELSE 
            DYU1(L,K)=0.     
          ENDIF
          IF(ISHDMF.EQ.2)THEN
            ! *** HMD WITH WALL EFFECTS
            IF(ICORDYU(L).EQ.2)THEN  
              DY2DZBR=1.+0.5*DYU(LS)/ZBRWALL  
              CSDRAG=0.16/((LOG(DY2DZBR))**2)  
              SLIPFAC=SLIPCO*SQRT(CSDRAG)  
              DYU1(L,K)=-2.*SLIPFAC*U(LS,K)/DYU(LS)  
            ENDIF  
            IF(ICORDYU(L).EQ.3)THEN  
              DY2DZBR=1.+0.5*DYU(L)/ZBRWALL  
              CSDRAG=0.16/((LOG(DY2DZBR))**2)  
              SLIPFAC=SLIPCO*SQRT(CSDRAG)  
              DYU1(L,K)=2.*SLIPFAC*U(L,K)/DYU(L)  
            ENDIF  
          ENDIF
          IF(ICORDXV(L).EQ.1)THEN
      ! *** DXV1 = dV/dX
           DXV1(L,K)=2.*SUB(L)*SUB(LW)*(V(L,K)-V(LW,K))/(DXV(L)+DXV(LW))
!           DXV1(L,K)=2.*(V(L,K)-V(LW,K))/(DXV(L)+DXV(LW))  !Hamrick formulation
          ELSE 
            DXV1(L,K)=0.     
          ENDIF
          IF(ISHDMF.EQ.2)THEN
            ! *** WALL EFFECTS
            IF(ICORDXV(L).EQ.2)THEN  
              DX2DZBR=1.+0.5*DXV(LW)/ZBRWALL  
              CSDRAG=0.16/((LOG(DX2DZBR))**2)  
              SLIPFAC=SLIPCO*SQRT(CSDRAG)  
              DXV1(L,K)=-2.*SLIPFAC*V(LW,K)/DXV(LW)  
            ENDIF  
            IF(ICORDXV(L).EQ.3)THEN  
              DX2DZBR=1.+0.5*DXV(L)/ZBRWALL  
              CSDRAG=0.16/((LOG(DX2DZBR))**2)  
              SLIPFAC=SLIPCO*SQRT(CSDRAG)  
              DXV1(L,K)=2.*SLIPFAC*V(L,K)/DXV(L)  
            END IF
          END IF
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
!$OMP END PARALLEL
C
C
      DO LL=1,NQWR
        DXU1(LIJ(IQWRU(LL),JQWRU(LL)),KQWRU(LL))=0.0
        DXV1(LIJ(IQWRU(LL),JQWRU(LL)),KQWRU(LL))=0.0
        DYU1(LIJ(IQWRU(LL),JQWRU(LL)),KQWRU(LL))=0.0
        DYV1(LIJ(IQWRU(LL),JQWRU(LL)),KQWRU(LL))=0.0
      ENDDO
      ! *** SXY = dU/dY + dV/dX     
      DO K=1,KC
        DO L=2,LA     
          SXY(L,K)=DYU1(L,K)+DXV1(L,K)     
        ENDDO  
      ENDDO    
C
C      DO K=1,KC      
C DO L=2,LA    
C   LS=LSC(L)  
C   LSW=LSWC(L)
C   TMPVAL=1.+SUB(L)+SVB(L)+SUB(LS)*SVB(L-1)   ! PMC     
C   DXYCOR=(DXYP(L)+SUB(L)*DXYP(L-1)+SVB(L)*DXYP(LS)     
C     & +SUB(LS)*SVB(L-1)*DXYP(LSW))/TMPVAL
C   DYU1(L,K)=DYU1(L,K)+DXV1(L,K)-0.5*     
C     & (SUB(L)*SUB(L-1)*(V(L,K)+V(L-1,K))*(DYV(L)-DYV(L-1))    
C     & +SVB(L)*SVB(LS )*(U(L,K)+U(LS ,K))*(DXU(L)-DXU(LS ))    
C     & )/DXYCOR      
C ENDDO 
C      ENDDO   
C
      IF(AHD.GT.0.0)THEN     
        ! *** CALCULATE SMAGORINSKY HORIZONTAL VISCOSITY 
!$OMP PARALLEL PRIVATE(TMPVAL,DSQR)
        DO K=1,KC     
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
         DO L=2,LA   
           TMPVAL=AHD*DXP(L)*DYP(L)
           DSQR=DXU1(L,K)*DXU1(L,K)+DYV1(L,K)*DYV1(L,K)+
     &         0.25*SXY(L,K)*SXY(L,K)
           AH(L,K)=AHO+TMPVAL*SQRT(DSQR)  
         ENDDO
        ENDDO  
!$OMP END PARALLEL

      ELSEIF(N.LT.10)THEN    
        ! *** ONLY NEED TO ASSIGN INITIALLY
        DO K=1,KC     
         DO L=2,LA   
           AH(L,K)=AHO      
         ENDDO
        ENDDO  
      ENDIF    
C
C *** DSLLC BEGIN BLOCK      
C **  CALCULATE HORIZONTAL DIFFUSION DUE TO WAVE BREAKING
C
      IF(ISWAVE.EQ.2)THEN    
        IF(WVLSH.GT.0.0.OR.WVLSX.GT.0.0)THEN      
          IF(N.LT.NTSWV)THEN 
            TMPVAL=FLOAT(N)/FLOAT(NTSWV)   
            WVFACT=0.5-0.5*COS(PI*TMPVAL)  
          ELSE 
            WVFACT=1.0
          ENDIF
          AHWVX=WVLSX*WVPRD*WVPRD   
!$OMP PARALLEL PRIVATE(DTMPH,DTMPX)          
          IF(ISDRY.GT.0)THEN
            DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
              DO L=2,LA
                IF(LMASKDRY(L))THEN  
                  DTMPH=WVDISP(L,K)**0.3333    
                  DTMPX=WVDISP(L,K)/HP(L)                    ! *** PMC HMP-->HP
                  AH(L,K)=AH(L,K)+WVFACT*(WVLSH*DTMPH*HP(L)
     &                     +AHWVX*DTMPX)         
                ENDIF
              ENDDO
!$OMP END DO NOWAIT
            ENDDO
          ELSE
            DO K=1,KC   
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
              DO L=2,LA 
                DTMPH=WVDISP(L,K)**0.3333    
                DTMPX=WVDISP(L,K)/HP(L)                      ! *** PMC HMP-->HP
                AH(L,K)=AH(L,K)+WVFACT*(WVLSH*DTMPH*HP(L)+AHWVX*DTMPX)  
              ENDDO
!$OMP END DO NOWAIT     
            ENDDO
          ENDIF
!$OMP END PARALLEL
        ENDIF  
      ENDIF    
C
C *** DSLLC END BLOCK 
C
      IF(N.EQ.2.AND.ISLOG.GT.0.AND.DEBUG)THEN     
      OPEN(1,FILE='AHDIFF.DIA')   
      CLOSE(1,STATUS='DELETE')    
      OPEN(1,FILE='AHDIFF.DIA')   
      DO L=2,LA     
        WRITE(1,1112)IL(L),JL(L),AH(L,KC)
      ENDDO  
      CLOSE(1)      
      ENDIF    
#ifdef key_mpi
      CALL COMMUNICATE_3D(AH)
#endif

C
C **  CALCULATE DIFFUSIVE MOMENTUM FLUXES  
C
!$OMP PARALLEL PRIVATE(LS,LN,LW,LE)
      DO K=1,KC
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=2,LA  
          LS=LSC(L)   
          LN=LNC(L)
          LW=LWEST(L)
          LE=LEAST(L)
          ! SANG'S CORRECTION
          FMDUX(L,K)=2.0*SUB(L)*
     &      (HP(L  )*AH(L  ,K)*DXU1(L  ,K)*DYP(L  )
     &      -HP(LW )*AH(LW ,K)*DXU1(LW ,K)*DYP(LW ))   

          FMDUY(L,K)=SVB(LN)*
     &      (DXU(LN)*HU(LN)*AH(LN,K)*SXY(LN,K)
     &      -DXU(L )*HU(L )*AH(L ,K)*SXY(L ,K))    

          FMDVY(L,K)=2.0*SVB(L)*
     &      (HP(L )*AH(L ,K)*DYV1(L ,K)*DXP(L )
     &      -HP(LS)*AH(LS,K)*DYV1(LS,K)*DXP(LS))

          FMDVX(L,K)=SUB(LE)*
     &      (DYV(LE )*HV(LE )*AH(LE ,K)*SXY(LE ,K)
     &      -DYV(L  )*HV(L  )*AH(L  ,K)*SXY(L  ,K))
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
!$OMP END PARALLEL


C
      ! *** TREAT THE NORTH & WEST WALL SLIPPAGE
      IF(ISHDMF.EQ.2)THEN
        DO L=2,LA
          LN=LNC(L)
          IF(SVBO(LN).LT.0.5)THEN
            DO K=1,KC
              DY2DZBR=1.+0.5*DYU(L)/ZBRWALL
              CSDRAG=0.16/((LOG(DY2DZBR))**2)
              SLIPFAC=SLIPCO*SQRT(CSDRAG)
              SXYLN=-2.*SLIPFAC*U(L,K)/DYU(L)
              FMDUY(L,K)=DXU(L)*HP(L)*AH(L,K)*(SXYLN-SXY(L ,K))
            ENDDO
          ENDIF

          IF(SUBO(LEAST(L)).LT.0.5)THEN
            DO K=1,KC
              DX2DZBR=1.+0.5*DXV(L)/ZBRWALL
              CSDRAG=0.16/((LOG(DX2DZBR))**2)
              SLIPFAC=SLIPCO*SQRT(CSDRAG)
              SXYEE=-2.*SLIPFAC*V(L,K)/DXV(L)
              FMDVX(L,K)=DYV(L)*HP(L)*AH(L,K)*(SXYEE-SXY(L,K))
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      ! *** ZERO BOUNDARY CELL MOMENTUM DIFFUSION
      DO LL=1,NBCS
        L=LBCS(LL)
        DO K=1,KC
          FMDUX(L,K)=0.0
          FMDUY(L,K)=0.0
          FMDVY(L,K)=0.0
          FMDVX(L,K)=0.0
        ENDDO
      ENDDO
      ! *** ZERO WITHDRAWAL/RETURN MOMENTUM DIFFUSION
      DO LL=1,NQWR
        FMDUX(LIJ(IQWRU(LL),JQWRU(LL)),KQWRU(LL))=0.0    
        FMDUY(LIJ(IQWRU(LL),JQWRU(LL)),KQWRU(LL))=0.0    
        FMDVX(LIJ(IQWRU(LL),JQWRU(LL)),KQWRU(LL))=0.0    
        FMDVY(LIJ(IQWRU(LL),JQWRU(LL)),KQWRU(LL))=0.0    
      ENDDO
C
      IF(N.EQ.2.AND.DEBUG)THEN  
        OPEN(1,FILE='AHD2.DIA')  
        CLOSE(1,STATUS='DELETE')  
        OPEN(1,FILE='AHD2.DIA')  
        DO L=2,LA  
          I=IL(L)  
          J=JL(L)  
          DO K=1,KC  
            WRITE(1,1111)N,I,J,K,FMDUX(L,K),FMDVY(L,K),FMDUY(L,K),  
     &          FMDVX(L,K),AH(L,K),DYU1(L,K),DYV1(L,K)  
          ENDDO  
        ENDDO  
        CLOSE(1)  
      ENDIF  
 1111 FORMAT(4I5,7E13.4)  
 1112 FORMAT(2I5,8E13.4)  
      RETURN  
      END  

