c
C****************************************************************************************
C****************************************************************************************
c****************************************************************************************
C
      SUBROUTINE PROJECTBLUE(U_TEMP,V_TEMP,DUBLUE,DVBLUE)
C
C **  SUBROUTINE PROJECTBLUE PROJECTS THE SURFACE ASSIMILATION OF BLUE INTO DEPTH USING EKAMN THEORY
C
C **  CREATED BY FEARGHAL O'DONNCHA ON 14 MARCH 2013
C----------------------------------------------------------------------------------------C
C
C
C
        USE GLOBAL
C
C*****************************************************************************************
C

        DOUBLE PRECISION  U_TEMP(LCM), V_TEMP(LCM),DUBLUE(LCM,KCM),DVBLUE(LCM,KCM), 
     &       DUS(LCM),DVS(LCM)
        REAL ZSIGMA,EKDEP(LCM),ALPHA

      IF(ISDYNSTP.EQ.0)THEN
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
      ELSE
        TIME=TIMESEC/TCTMSR
      ENDIF
      
      DO K =1,KC
         DO L = 2,LA
            DUBLUE(L,K) = 0.
            DVBLUE(L,K) = 0.
         END DO
       END DO   


        DO K =1,KS
           ZSIGMA = (K/(KC*1.))-1.
           DO L =2,LA
           ekdep(L) = SQRT(  (2.*AV(L,KS))  /CF )

! ekdep = EKMAN DEPTH = SQRT(2AV/f)   
!     AV(L,K) = model computed viscosity
!     cf = constant coriolis parameter(1/s) =2*7.29E-5*(SIN(LAT)); galway bay lat = 53.2

           ALPHA = (ZSIGMA*HU(L))/EKDEP(L)        
!           dus = U(L,KC) - U_TEMP(L)        ! surface shear, i.e. U(model) - U(fn(model,sensor))
!           dvs = V(L,KC) - V_TEMP(L)

           dus(L) = U_TEMP(L) - U(L,KC)        ! surface shear, i.e. U(model) - U(fn(model,sensor))
           dvs(L) = V_TEMP(L) - V(L,KC)
          DUBLUE(L,K) =  EXP(ALPHA) * ( (dus(L)*COS(ALPHA)) - (dvs(L)*SIN(ALPHA)))
          DVBLUE(L,K) =  EXP(ALPHA) * ( (dus(L)*SIN(ALPHA)) + (dvs(L)*COS(ALPHA)))
            END DO
        END DO



! write computed EKMAN profile to file at every assimilation point

        DO MLTM = 1,MLTMSR
           OPEN(123,FILE='DUBLUE'//CNTMSR(MLTM)//'.csv',status='unknown',position='append')
           I = ILTMSR(MLTM)
           J = JLTMSR(MLTM)
           L = LIJ(I,J)
           write(123,123)I,J,TIME,(DUBLUE(L,K),K=1,KC),EKDEP(L),DUS(L),DT,DT2
           CLOSE(123)

           OPEN(124,FILE='DVBLUE'//CNTMSR(MLTM)//'.csv',status='unknown',position='append')
      
           write(124,123)I,J,TIME,(DVBLUE(L,K),K=1,KC),EKDEP(L),DVS(L)
           CLOSE(124)
           END DO
 123       FORMAT(2I5,F12.5,52E12.4,2f6.1)

        RETURN
        END
