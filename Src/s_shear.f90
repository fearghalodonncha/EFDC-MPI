SUBROUTINE SEDZLJ_SHEAR
  
  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER::L
  INTEGER::M1,M2
  INTEGER::FZONE
  !PT: All real values are explicitly written in DOUBLE PRECISION 7/16/08.
  DOUBLE PRECISION::MMW,SIGMAWV,JJW
  DOUBLE PRECISION::VELMAG,VELANG,DELW,APROUGH
  DOUBLE PRECISION::UTMP,VTMP
  DOUBLE PRECISION::WVLENGTH,WVANGLE,WFTIM,EXCURSION
  DOUBLE PRECISION::FC1,FC2,FWINDSQ
  DOUBLE PRECISION::FWINDS,FWINDD
  DOUBLE PRECISION::TDIFF,WTM1,WTM2,AVGDEPTH
  DOUBLE PRECISION,DIMENSION(LCM)::ZBTEMP
  REAL:: FC3  
  ! CALCULATES Wave and Current Shear Stress Based on Log-Law of Cristofferson Jonsson
  ! 
  ! REVISION DATE :  May 24, 2007
  ! Craig Jones and Scott James
  !**************************************************************************
  ! Check to see if we've set a constant Tau
  
  IF(TAUCONST==0)THEN
     !
     !******************************************************
     IF(ISWNWAVE==1)THEN
        !
        !**************************************************************************
        ! Wind Wave Fetch
        !
        ! Convert wind input into current wind info for wind-driven wave calcs
        IF(ISDYNSTP==0)THEN  
           WFTIM=DT*FLOAT(N)/TCWSER(1)+TBEGIN*(TCON/TCWSER(1))  
        ELSE  
           WFTIM=TIMESEC/TCWSER(1)  
        ENDIF
        M1=MWTLAST(1)  
        M2=M1+1  
        DO WHILE(TWSER(M2,1)<WFTIM)  
           M1=M2
           M2=M1+1  
        ENDDO
        MWTLAST(1)=M1 
        TDIFF=TWSER(M2,1)-TWSER(M1,1)  
        WTM1=(TWSER(M2,1)-WFTIM)/TDIFF  
        WTM2=(WFTIM-TWSER(M1,1))/TDIFF 
        FWINDS=WTM1*WINDS(M1,1)+WTM2*WINDS(M2,1)
        IF(FWINDS>1.0)THEN
           IF(ABS(WINDD(M1,1)-WINDD(M2,1))<180.0)THEN
              FWINDD=WTM1*WINDD(M1,1)+WTM2*WINDD(M2,1)
           ELSE
              IF(WINDD(M1,1).GT.WINDD(M2,1))THEN
                 FWINDD=WTM1*WINDD(M1,1)+WTM2*(WINDD(M2,1)+360.0)
              ELSE
                 FWINDD=WTM1*(WINDD(M1,1)+360)+WTM2*WINDD(M2,1)
              ENDIF
              IF(FWINDD>=360.0)FWINDD=FWINDD-360.0 
           ENDIF
           ! Convert wind into direction it is blowing "from"
           IF(FWINDD<=180.0)THEN  
              FWINDD=FWINDD+180.0  
              IF(FWINDD==360.0)FWINDD=0.0         
           ELSE  
              FWINDD=FWINDD-180.0  
              IF(FWINDD==360.0)FWINDD=0.0 
           ENDIF
           ! Calculate which of the 8 wind zones (FWZONE) the wind is coming from
           ! Also the Waveangle CCW from East
           IF(FWINDD>=337.5.OR.FWINDD<22.5)THEN
              FZONE=1
              WVANGLE=4.712
           ELSEIF(FWINDD>=22.5.AND.FWINDD<67.5)THEN
              FZONE=2
              WVANGLE=3.927
           ELSEIF(FWINDD>=67.5.AND.FWINDD<112.5)THEN
              FZONE=3
              WVANGLE=3.142
           ELSEIF(FWINDD>=112.5.AND.FWINDD<157.5)THEN
              FZONE=4
              WVANGLE=2.356
           ELSEIF(FWINDD>=157.5.AND.FWINDD<202.5)THEN
              FZONE=5
              WVANGLE=1.571
           ELSEIF(FWINDD>=202.5.AND.FWINDD<247.5)THEN
              FZONE=6
              WVANGLE=0.7854
           ELSEIF(FWINDD>=247.5.AND.FWINDD<292.5)THEN
              FZONE=7
              WVANGLE=0.
           ELSEIF(FWINDD>=292.5.AND.FWINDD<337.5)THEN
              FZONE=8
              WVANGLE=5.4978
           ENDIF
           !Calculate Domain Average Depth
           !Needs to be calculated along each fetch
           !This is sufficient for small systems
           AVGDEPTH=SUM(HP(2:LA))/FLOAT(LA-1)
           FWINDSQ=FWINDS*FWINDS
        ENDIF
           
           !Calculate wave height, period, orbital velocity, and length
           DO L=2,LA
              FC1=(FWINDSQ/9.8)*0.283*TANH(0.530*(9.8*AVGDEPTH/FWINDSQ)**0.75)
              FC2=TANH(0.0125*(9.8*FWDIR(L,FZONE)/FWINDSQ)**0.42/TANH(0.530*(9.8*AVGDEPTH/FWINDSQ)**0.75))   
              FC3=FC1*FC2  
              FWVHT(L)=MIN(HP(L),FC3) !taking min of real/doubleprecision causes issues in AIX; 
              !FWVHT(L)=MIN(0.6*HP(L),FWVHT(L))
              
              FC1=(FWINDS/9.8)*7.54*TANH(0.833*(9.8*AVGDEPTH/FWINDSQ)**0.375)
              FC2=TANH(0.077*(9.8*FWDIR(L,FZONE)/FWINDSQ)**0.25/TANH(0.833*(9.8*AVGDEPTH/FWINDSQ)**0.375))   
              FWVTP(L)=FC1*FC2
              WVFREQ(L)=2.0*PI/FWVTP(L)
              
              !FC1=(2.0*PI/FWVTP(L))**2*HP(L)/9.8
              !FC2=FC1+1.0/(1.0+0.6522*(FC1)+0.4622*(FC1)**2+0.0864*(FC1)**4+0.0675*(FC1)**5)
              WVLENGTH=FWVTP(L)*SQRT(9.8*HP(L)/FC2)
              WVLENGTH=MAX(1.d0,WVLENGTH)
              !WVORBIT(L)=PI*FWVHT(L)/(FWVTP(L)*SINH(HP(L)*2.0*PI/WVLENGTH))
              WVORBIT(L)=MAX(0.01d0,PI*FWVHT(L)/(FWVTP(L)*SINH(HP(L)*2.0*PI/WVLENGTH)))
           ENDDO
           !*************************************************************************
           !Read in EFDC STWAVE Wave Field
           
      ELSEIF(ISWNWAVE==2)THEN
           
           NWVCOUNT=NWVCOUNT+DT/3600
           IF(NWVCOUNT==STWVTIM)THEN
              NWVCOUNT=0
              STINC=STINC+1
              DO L=2,LA
                 IF(STINC>STWVNUM)EXIT
                 IF(STWVTP(L,STINC)>0.0.AND.ISCDRY(L).EQ.0)THEN
                    WVFREQ(L)=2.0*PI/STWVTP(L,STINC)
                    FC3 = STWVHT(L,STINC)
                    FWVHT(L)=MIN(HP(L),FC3)
!                    FC1=(2.0*PI/STWVTP(L,STINC))**2*HP(L)/9.8
!                    FC2=FC1+1.0/(1.0+0.6522*(FC1)+0.4622*(FC1)**2+0.0864*(FC1)**4+0.0675*(FC1)**5)
!                    WVLENGTH=STWVTP(L,STINC)*SQRT(9.8*HP(L)/FC2)
                    FC1=9.8*STWVTP(L,STINC)**2
                    FC2=2.*PI*SQRT(TANH(4.*PI**2*HP(L)/(STWVTP(L,STINC)**2*9.8)))
                    WVLENGTH=FC1/FC2
                    EXCURSION=FWVHT(L)/(2.*SINH((2.*PI/WVLENGTH)*HP(L)))
                    WVORBIT(L)=EXCURSION*WVFREQ(L)
                    WVORBIT(L)=MAX(0.01d0,WVORBIT(L))
                    WVANG(L)=STWVDR(L,STINC)
                 ELSE
                    WVFREQ(L)=0.0
                    WVORBIT(L)=0.0
                    WVANG(L)=0.0
                 ENDIF
              ENDDO
           ENDIF
       ENDIF    
     !*************************************************************************
     !Begin shear stress calculations
     !Set up roughness, velocities, and angles
     
     DO L=2,LA
        ZBTEMP(L)=ZBR(L)
        IF(ZBSKIN.EQ.0.AND.NSEDFLUME.GT.0)THEN
           IF(D50AVG(L).LT.D50(1))THEN
              ZBTEMP(L)=D50(1)/1e6
           ELSE
              ZBTEMP(L)=D50AVG(L)/1e6
           ENDIF
        ELSEIF(ZBSKIN.GT.0.AND.NSEDFLUME.GT.0)THEN
           ZBTEMP(L)=ZBSKIN/1e6
        ENDIF
     ENDDO
     
     DO L=2,LA
        IF(LMASKDRY(L).AND.HP(L).GT.0.5)THEN
          IF(ISWAVE.EQ.2) ZBTEMP(L)=0.0001
          KN(L)=30.0*ZBTEMP(L)
          ! Calculate Average Velocity Magnitude in cm/s
          UTMP=100.0*STCUV(L)*(UHE(LEAST(L))+UHE(L))/(HU(LEAST(L))+HU(L))+1.0E-12
          VTMP=100.0*STCUV(L)*(VHE(LNC(L))+VHE(L))/(HV(LNC(L))+HV(L))
          VELMAG=SQRT(UTMP**2+VTMP**2)
          ! Calculate Initial Friction Factors
          FC(L)=(0.42/LOG(HP(L)/(2.0*ZBTEMP(L))))**2	
          ! Current only Shear Stress
          IF(ISWNWAVE==0.AND.UWVSQ(L)==0.0.OR.ISWAVE==0)THEN
             TAU(L)=FC(L)*VELMAG**2
          ELSE
             !Calculate Combined Wave Friction Factor
             ! Calculate Current Angle CCW From X axis
             IF(UTMP>0.0.AND.VTMP>0.0)THEN
                VELANG=ATAN(VTMP/UTMP)
             ELSEIF(UTMP<0.0)THEN
                VELANG=ATAN(VTMP/UTMP)+PI
             ELSEIF(UTMP>0.0.AND.VTMP<0.0)THEN
                VELANG=2*PI+ATAN(VTMP/UTMP)
             ELSEIF(UTMP==0.0)THEN
                VELANG=SIGN(0.5d0*PI,VTMP)
             ENDIF
             ! Set Orbital velocity in m/s and waveangle and frequency
             IF(ISWAVE.GT.0)THEN
                 WVFREQ(L)=WVFRQ
                 WVORBIT(L)=SQRT(UWVSQ(L))
                 WVANG(L)=WACCWE(L)
             ELSEIF(ISWNWAVE==1)THEN
                WVANG(L)=WVANGLE
             ENDIF
             ! Calculate wave friction factor
             FWW(L)=2.0*(0.0747*(KN(L)*WVFREQ(L)/WVORBIT(L)))**0.66
             SIGMAWV=FC(L)/FWW(L)*(VELMAG/(WVORBIT(L)*100.0))**2
             MMW=SQRT(1.0+SIGMAWV**2+2.0*SIGMAWV*ABS(COS(VELANG-WVANG(L))))
             JJW=WVORBIT(L)/(KN(L)*WVFREQ(L))*SQRT(MMW*FWW(L)/2.0)
             FWW(L)=MMW*0.15/JJW
             !Calculate wave boundary layer info
             DELW=KN(L)*0.273*SQRT(JJW)
             APROUGH=30.0*DELW*EXP(-5.62*DELW/KN(L)*SQRT(SIGMAWV/MMW))
             !Calculate new current friction factor
             FC(L)=2.0*(1.0/(2.38*LOG(30.0*HP(L)/(2.718*KN(L)))-2.38*LOG(APROUGH/KN(L))))**2
             !Iterate once more
             SIGMAWV=FC(L)/FWW(L)*(VELMAG/(WVORBIT(L)*100.0))**2
             MMW=SQRT(1.0+SIGMAWV**2+2.0*SIGMAWV*ABS(COS(VELANG-WVANG(L))))
             JJW=WVORBIT(L)/(KN(L)*WVFREQ(L))*SQRT(MMW*FWW(L)/2.0)
             FWW(L)=MMW*0.15/JJW
             !Calculate total wave and current shear stress (dynes/cm^2)
             TAU(L)=0.5*FWW(L)*WVORBIT(L)**2*10000.0*MMW
!             TAUB(L)=0.1*TAU(L) !PT: conversion from dynes/cm^2 to Pa Don't overwrite EFDC's total bed shear stress
             USTAR(L)=SQRT(TAUB(L)/1000.0) !only true because USTAR=SQRT(Tau/RhoH2O) and RhoH2O is 1000 kg/cm^3. 
          ENDIF
        ENDIF
     ENDDO
     !*****************************
     
  ELSE !Set constant tau is TAUCONST (dynes/cm^2) is greater than 0
     DO L=2,LA
        TAU(L)=TAUCONST
 !       TAUB(L)=0.1*TAU(L) !Don't overwrtie EFDC's shear stress.
        USTAR(L)=SQRT(TAUB(L)/1000.)
     ENDDO
  ENDIF
  !************************
  RETURN  
END SUBROUTINE SEDZLJ_SHEAR
