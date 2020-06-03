      SUBROUTINE CALAVB (ISTL_)  
C  
C **  SUBROUTINE CALAV CALCULATES VERTICAL VISCOSITY AND DIFFUSIVITY  
C **  USING GLAPERIN ET AL'S MODIFICATION OF THE MELLOR-YAMADA MODEL  
C **  (NOTE AV, AB, AND AQ ARE ACTUALLY DIVIDED BY H)  
C **  IF ISGA=1 VALUES ARE GEOMETRIC AVERAGES WITH THE PREVIOUS VALUES  
C CHANGE RECORD  
C  ADDED DRYCELL BYPASS AND CONSISTENT INITIALIZATION OF DRY VALUES  
C  
      USE GLOBAL  
      USE OMP_LIB
	IMPLICIT NONE
	INTEGER::L,K,LS,ISTL_
	REAL::QQIMAX,RIQMIN,RIQMAX,RIQ,SFAV,SFAB,ABTMP,AVTMP
      REAL::BVF,SH2,SHEARVEL,SURFTAU,SURFVEL
!     New variable needed to find the the depth to the thermocline
      REAL::ZTHERM(LCM),ZMET(LCM),DENMAX(LCM)
      INTEGER::KZTHERM(LCM),KZMET(LCM)
      REAL::P2,alpha,Knot,Knottherm,Knotmet,BVFZ,Kv0  !SCJ deleted CP1,CP2,CP3,
      REAL::RIQZtherm,RIQZmet,RIQZ,BVFZtherm,BVFZmet,KZTH,KZM
      
      REAL::RIQprint

C  
C   SHTOP    =      0.4939  
C   SHBOT    =     34.6764  
C   SMTOP1   =      0.3933  
C   SMTOP2   =      7.8464  
C   SMBOT1   =     34.6764  
C   SMBOT2   =      6.1272  
C   RLIMIT   =      0.0233  
C   SHMIN    =      0.0934  
C   SMMIN    =      0.1099  
C   SHMAX    =      5.2073  
C   SMMAX    =      4.9639  
C  
      QQIMAX=1./QQMIN  
      AVMAX=AVO  
      ABMAX=ABO  
      AVMIN=10.  
      ABMIN=10.  
      IF(ISFAVB==3.OR.ISFAVB==4)THEN !Implement Lozovatsky Ri range
        RIQMIN=0.0 !was 0.25
        RIQMAX=2.0
      ELSE  
        RIQMIN=-0.023  
        RIQMAX=0.28  
      ENDIF
      DO K=1,KC  
        DO L=1,LC  
          IF(IMASKDRY(L).EQ.1)THEN  
            AV(L,K)=AVO*HPI(L)  
            AB(L,K)=ABO*HPI(L)  
          ENDIF  
        ENDDO  
      ENDDO  
      IF(ISFAVB.EQ.0)THEN  
!$OMP PARALLEL PRIVATE(RIQ,SFAV,SFAB,AVMAX,ABMAX,AVMIN,ABMIN)
        DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              QQI(L)=1./QQ(L,K)  
              QQI(L)=MIN(QQI(L),QQIMAX)  

              RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)  
     &            *(B(L,K+1)-B(L,K))*QQI(L)  
              RIQ=MAX(RIQ,RIQMIN)  
              RIQ=MIN(RIQ,RIQMAX)  
C  
C      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))  
C      SFAB=0.5/(1.+36.*RIQ)  
C      SFAV=0.3933*(1.+7.8464*RIQ)/((1.+34.6764*RIQ)*(1.+6.1272*RIQ))  
C      SFAB=0.4939/(1.+34.6764*RIQ)  
C SFAV AND SFAB ABOVE REPLACE BY KANATHA AND CLAYSON  
C  
              SFAV=0.3920*(1.+8.6736*RIQ)/((1.+30.192*RIQ)*(1.+
     &            6.1272*RIQ))  
              SFAB=0.4939/(1.+30.192*RIQ)  
             AB(L,K)=AVCON*SFAB*DML(L,K)*HP(L)*QQSQR(L,K)+ABO  
             AV(L,K)=AVCON*SFAV*DML(L,K)*HP(L)*QQSQR(L,K)+AVO  
              AVMAX=MAX(AVMAX,AV(L,K))  
              ABMAX=MAX(ABMAX,AB(L,K))  
              AVMIN=MIN(AVMIN,AV(L,K))  
              ABMIN=MIN(ABMIN,AB(L,K))  
              AV(L,K)=AV(L,K)*HPI(L)  
              AB(L,K)=SCB(L)*AB(L,K)*HPI(L)  
            ENDIF  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
!$OMP END PARALLEL
      ELSEIF (ISFAVB.EQ.1)THEN 
!$OMP PARALLEL PRIVATE(RIQ,SFAV,SFAB,ABTMP,AVTMP,AVMAX,ABMAX,AVMIN,ABMIN)
        DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              QQI(L)=1./QQ(L,K)  
              QQI(L)=MIN(QQI(L),QQIMAX)  
              RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)  
     &            *(B(L,K+1)-B(L,K))*QQI(L)  
              RIQ=MAX(RIQ,RIQMIN)  
              RIQ=MIN(RIQ,RIQMAX)  
C  
C      SFAV=0.4*(1.+8.*RIQ)/((1.+36.*RIQ)*(1.+6.*RIQ))  
C      SFAB=0.5/(1.+36.*RIQ)  
C      SFAV=0.3933*(1.+7.8464*RIQ)/((1.+34.6764*RIQ)*(1.+6.1272*RIQ))  
C      SFAB=0.4939/(1.+34.6764*RIQ)  
C SFAV AND SFAB ABOVE REPLACE BY KANATHA AND CLAYSON  
C  
              SFAV=0.3920*(1.+8.6736*RIQ)/((1.+30.192*RIQ)*(1.+
     &            6.1272*RIQ))  
              SFAB=0.4939/(1.+30.192*RIQ)  
              ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*QQSQR(L,K)+ABO  
              AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*QQSQR(L,K)+AVO  
              AVMAX=MAX(AVMAX,AVTMP)  
              ABMAX=MAX(ABMAX,ABTMP)  
              AVMIN=MIN(AVMIN,AVTMP)  
              ABMIN=MIN(ABMIN,ABTMP)  
              AV(L,K)=0.5*(AV(L,K)+AVTMP*HPI(L))  
              AB(L,K)=SCB(L)*0.5*(AB(L,K)+ABTMP*HPI(L))  
            ENDIF  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
!$OMP END PARALLEL

      ELSEIF(ISFAVB.EQ.2)THEN  
!$OMP PARALLEL PRIVATE(RIQ,SFAV,SFAB,ABTMP,AVTMP,AVMAX,ABMAX,AVMIN,ABMIN)
        DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              QQI(L)=1./QQ(L,K)  
              QQI(L)=MIN(QQI(L),QQIMAX)  
              RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)  
     &            *(B(L,K+1)-B(L,K))*QQI(L)  
              RIQ=MAX(RIQ,RIQMIN)  
              RIQ=MIN(RIQ,RIQMAX)  
C  
C SFAV AND SFAB ABOVE REPLACE BY KANATHA AND CLAYSON  
C  
              SFAV=0.3920*(1.+8.6736*RIQ)/((1.+30.192*RIQ)*(1.+
     &            6.1272*RIQ))  
              SFAB=0.4939/(1.+30.192*RIQ)  
              ABTMP=AVCON*SFAB*DML(L,K)*HP(L)*QQSQR(L,K)+ABO  
              AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*QQSQR(L,K)+AVO  
              AVMAX=MAX(AVMAX,AVTMP)  
              ABMAX=MAX(ABMAX,ABTMP)  
              AVMIN=MIN(AVMIN,AVTMP)  
              ABMIN=MIN(ABMIN,ABTMP)  
              AV(L,K)=SQRT(AV(L,K)*AVTMP*HPI(L))  
              AB(L,K)=SCB(L)*SQRT(AB(L,K)*ABTMP*HPI(L))  
            ENDIF  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
!$OMP END PARALLEL
      ELSEIF(ISFAVB.EQ.3)THEN  !R Arafin formulation 1
        DO K=1,KS  !KS=KC-1 second layer down
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              QQI(L)=1./QQ(L,K)  
              QQI(L)=MIN(QQI(L),QQIMAX)  
            ENDIF  
          ENDDO  
          DO L=2,LA  
            IF(LMASKDRY(L))THEN
!Brunt-Vaisala Frequency: N^2 =-g/rho0 * (del rho)/(del z)
!(del rho)/rho0 = B(L,K+1)-B(L,K)
!1/del z = HPI(L)*DZIG(K)
              BVF=-GP*(B(L,K+1)-B(L,K))*HPI(L)*DZIG(K)
!Vertical shear: [(del u)/(del z)]^2 + [(del v)/(del z)]^2
              SH2=((U(L,K+1)-U(L,K))*HPI(L)*DZIG(K))**2+
     &            ((V(L,K+1)-V(L,K))*HPI(L)*DZIG(K))**2
              RIQ=BVF/MIN(SH2+1.0E-6,RIQMAX) !Richardson number is Brunt-Vaisala frequency divided by vertical shear  
              RIQ=MAX(RIQ,RIQMIN)  !Implement floor of Ri
              RIQ=MIN(RIQ,RIQMAX)  !Implement ceiling of Ri
!Shear velocity calculation: u*=sqrt(tau/rho)=sqrt(mu*sqrt(SH2)/rho)
!Water viscosity:            mu=1.003E-3 kg/m-s
!Water density:              rho=1000.0*(B(L,K)+1.0)
              SHEARVEL=
     &SQRT(1.003E-3*SQRT(SH2)/(1000.0*(B(L,K)+1.0)))
!Vertical diffusivity: von Karman constant is VKC
!                      Depth, z, is HP(L)*0.5*(Z(K+1)+Z(K))
!                      p is (2.0/3.0)
!                      Ricr is 0.1
!                      Ribeta is 0.05
!                      r is 1.0
              IF(HP(L)*(1.0-ZZ(K))<30.0)THEN !Lozovatsky formulation for depths less that 30 m
                SFAB=VKC*SHEARVEL*HP(L)*(1.0-ZZ(K))/
     &             (1.0+RIQ/0.1)**(2.0/3.0)/(1.0+RIQ/0.05)**1.0
                RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)  !Old Mellor-Yamada formulation for SFAV
     &              *(B(L,K+1)-B(L,K))*QQI(L)  
                RIQ=MAX(RIQ,RIQMIN)  
                RIQ=MIN(RIQ,RIQMAX)
                SFAV=0.3920*(1.+8.6736*RIQ)/
     &               ((1.+30.192*RIQ)*(1.+6.1272*RIQ))
              ELSE !Mellor-Yamada formulation for depths greater than 30 m
                RIQ=-GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(K)  
     &            *(B(L,K+1)-B(L,K))*QQI(L)  
                RIQ=MAX(RIQ,RIQMIN)  
                RIQ=MIN(RIQ,RIQMAX)  
                SFAV=0.3920*(1.+8.6736*RIQ)/((1.+30.192*RIQ)*(1.+
     &            6.1272*RIQ))
                SFAB=0.4939/(1.+30.192*RIQ) 
              ENDIF

              IF(K==KS)THEN !Special logic for top model layer not otherwise calculated. In the Mellor-Yamada version, it is just ABO*HPI(L)
                  SURFTAU=SQRT(TSX(L)**2 + TSY(L)**2)           !RRA: TSX(L) and TSY(L) computed in CALTSXY
                  SURFVEL=SQRT(SURFTAU/(1000.0*(B(L,K+1)+1.0))) !RRA: SURFVEL is frictional velocity U* at surface
                  ABTMP=AVCON*VKC*SURFVEL*HP(L)*(1.0-ZZ(K+1))/
     &                  (1.0+RIQ/0.1)**(2.0/3.0)/(1.0+RIQ/0.05)**1.0+ABO
                  AB(L,K+1)=SCB(L)*0.5*(AB(L,K+1)+ABTMP*HPI(L))
              ENDIF                  
              ABTMP=AVCON*SFAB+ABO  
              AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*QQSQR(L,K)+AVO  
              AVMAX=MAX(AVMAX,AVTMP)  
              ABMAX=MAX(ABMAX,ABTMP)  
              AVMIN=MIN(AVMIN,AVTMP)  
              ABMIN=MIN(ABMIN,ABTMP)  
              AV(L,K)=0.5*(AV(L,K)+AVTMP*HPI(L))  !Harmonic (time average) time smoothing
              AB(L,K)=SCB(L)*0.5*(AB(L,K)+ABTMP*HPI(L))  !Harmonic (time averge) time smoothing
!      IF(MOD(L,1000)==0)write(6,'(3(i5,1x),6(e9.2,1x))') !Debugging
!     &N,L,K,
!     &AB(L+1,K),AB(L,K),
!     &AB(L,K+1),RIQ
            ENDIF  
          ENDDO  
        ENDDO
      ELSEIF(ISFAVB.EQ.4)THEN
!Variation of density B(L,K+1)-B(L,K) 
!Variation of vertical velocities U(L,K+1)-U(L,K) and V(L,K+1)-V(L,K)
!Change in sigma depth (normalized Z) is DZC(K)
!Change in  depth Z (dimensional) is HP(L)*DZC(K)
!RRA: Need to define depth, z correctly 
!Ztherm is the depth of maximum gradient
!Zmeta is the boundary depth between metalimnion and hypolimnion at which the density gradient equals to threshold value 10-5 kg m-4
        IF (PARTID == MASTER_TASK) THEN
          OPEN(9876,FILE="ZTHERM.OUT",POSITION='APPEND')
          IF(MOD(N-1,86400)==0)THEN
            WRITE(9876,'("Time step: ",I8)')N
            WRITE(9876,
     &  '("  I   J Ktherm HP   Ztherm  Zmet Kmet  AVKT  ABKT ")')
          ENDIF
        END IF
        DENMAX(2:LA)=0.0 !initialize maximum density change to zero
        KZTHERM(2:LA)=KS !initialize the thermocline to start at layer KS (one below top layer)
        ZTHERM(2:LA)=HP(2:LA)*(1.0-ZZ(KS)) !initialize thermocline depth to the depth of layer KS
        KZMET(2:LA)=2 !initialize the metalimnion to begin at layer 2
        ZMET(2:LA)=HP(2:LA)*(1.0-ZZ(2)) !initialize the metalimnion depth to the depth of layer 2
        DO L=2,LA
          DO K=KS,2,-1
            IF(DENMAX(L)<B(L,K)-B(L,K+1))THEN !find the depth of the maximum density gradient (assumes stable stratification)
              DENMAX(L)=B(L,K)-B(L,K+1)
              ZTHERM(L)=HP(L)*(1.0-ZZ(K))
              KZTHERM(L)=K !Note this KZTHERM is the layer (K) containing the thermocline and not K(Ztherm), which is calculated below as KZTH
            ENDIF
          ENDDO
          IF(HP(L)>30.0)THEN !deep water
            IF(ZTHERM(L)<4.0)THEN !limit thermocline to deeper than 4 m
              ZTHERM(L)=4.0
              DO K=KC,1,-1
                IF(HP(L)*(1.0-ZZ(K))>4.0)THEN
                  KZTHERM(L)=MAX(K,2)
                  EXIT
                ENDIF
              ENDDO
            ELSEIF(ZTHERM(L)>30.0)THEN !limit thermocline to shallower than 30 m
              ZTHERM(L)=30.0
              DO K=1,KC
                IF(HP(L)*(1.0-ZZ(K))<30.0)THEN
                  KZTHERM(L)=MAX(K,2)
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ELSEIF(HP(L)>12.0)THEN !shallow water
            IF(ZTHERM(L)<4.0)THEN !limit thermocline to deeper than 4 m
              ZTHERM(L)=4.0
              DO K=KC,1,-1
                IF(HP(L)*(1.0-ZZ(K))>4.0)THEN
                  KZTHERM(L)=MAX(K,2)
                  EXIT
                ENDIF
              ENDDO
            ELSEIF(ZTHERM(L)>12.0)THEN !limit thermocline to shallower than 12 m
              ZTHERM(L)=12.0
              DO K=1,KC
                IF(HP(L)*(1.0-ZZ(K))<12.0)THEN
                  KZTHERM(L)=MAX(K,2)
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ELSE
            CONTINUE !Do not limit thermocline depths if the water depth is less than 12 m
          ENDIF  
          DO K=KZTHERM(L)-1,2,-1 !calculate the metalimnion depth.
            IF(1000.0*(B(L,K)-B(L,K-1))/(DZC(K)*HP(L))<1.0E-5)THEN !Is the density gradient < 1E-5 kg/m^4
              ZMET(L)=HP(L)*(1.0-ZZ(K))
              KZMET(L)=K !Note this KZMET is the layer (K) containing the metalimnion and not K(Zmet), which is calculated below as KZM
             EXIT
            ENDIF
          ENDDO
        ENDDO       
        DO K=1,KS  !KS=KC-1 (second layer from the top)
          DO L=2,LA  
            IF(LMASKDRY(L))THEN  
              QQI(L)=1./QQ(L,K)  
              QQI(L)=MIN(QQI(L),QQIMAX)  
            ENDIF  
          ENDDO
        ENDDO
        DO L=2,LA  
          IF(LMASKDRY(L))THEN  
!Surface frictional velocity computation: 2 approaches, either default from EFDC or from Vincon-Leite.
!Water viscosity:           mu      = 1.003E-3 kg/m-s
!Water density:             rho     = 1000.0*(B(L,K)+1.0)
!Air density: 			  rho_air = 1.29*(B(L,K)+1.0)
!Wind drag coefficient, CD: CD      = 1.8 E-3      

!Surface frictional velocity computation as default in EFDC.
!              SURFTAU=SQRT(TSX(L)**2 + TSY(L)**2)           !RRA: TSX(L) and TSY(L) computed in CALTSXY
!              SURFVEL=SQRT(SURFTAU/(1000.0*(B(L,K+1)+1.0))) !RRA: SURFVEL is frictional velocity U* at surface
!Surface frictional velocity computation following Vincon-Leite. 
            SURFVEL=
     &SQRT(1.8E-3*(1.29*(B(L,KS)+1.0))/(1000.0*(B(L,KC)+1.0)))*WINDST(L) !RRA: WINDST(L)is computed in CALTSXY This is calculated only at the surface layer
            DO K=KS,1,-1
!Richardson number, RIQ computation as default in EFDC (Brunt-Vaisala frequency divided by vertical shear)
!Vertical shear: [(del u)/(del z)]^2 + [(del v)/(del z)]^2
!              SH2=((U(L,K+1)-U(L,K))*HPI(L)*DZIG(K))**2+
!     &            ((V(L,K+1)-V(L,K))*HPI(L)*DZIG(K))**2
!              RIQ=BVF/MIN(SH2+1.0E-6,RIQMAX) !Richardson number is Brunt-Vaisala frequency divided by vertical shear
!Vincon-Leite Richardson number, RIQ, is independent of shear velocity
! Calibration parameter CP1: gamma
              CP1 = 570.0 !seconds
! Calibration parameter CP2: delta
              CP2 = 0.015
! Calibration parameter CP3: sigma
              CP3 = 9.46
! Calibration parameter alpha:
              alpha = 0.173
! Calibration parameter P2:
              P2 = 0.605
!Brunt-Vaisala Frequency: N^2 =-g/rho0 * (del rho)/(del z)
!(del rho)/rho0 = B(L,K+1)-B(L,K)
!1/del z = HPI(L)*DZIG(K)
              BVF=GP/(1.0+B(L,K))*ABS(B(L,K+1)-B(L,K))*HPI(L)*DZIG(K) !Brunt-Vaisala frequency squared at depth Z (layer K): BVF(Z) = g/rho*|delta rho|/|delta z|
              
!RRA: Shear velocity (Internal) computation is redundant because Richardson number is independent of shear velocity according to Vincon-Leite.
!u*=sqrt(tau/rho)=sqrt(mu*sqrt(SH2)/rho)
!Water viscosity:            mu=1.003E-3 kg/m-s
!Water density:              rho=1000.0*(B(L,K)+1.0)
!SHEARVEL=SQRT(1.003E-3*SQRT(SH2)/(1000.0*(B(L,K)+1.0)))
!Eddy Diffusivity K0 in neutral condition
!RRA: Note RIQZ (Richardson number) at metalimnion and hypolimnion levels required in K(Ztherm) and K(Zmet) computations
!BVFZ (Brunt-Vaisala frequency) required at Ztherm and Zmet depths
             IF(HP(L)*(1.0-ZZ(K))<=ZTHERM(L))THEN !From Surface to Thermocline
               RIQ=
     &BVF*(CP1*CP2)**2*(EXP(2.*HP(L)*(1.0-ZZ(K))/CP1/SURFVEL)) !Ri(Z)
!RRA: Vincon-Leite didn't mention any upper or lower limit for Richardson number, so may be floor/ceiling of Ri not needed.
               RIQ=MAX(RIQ,RIQMIN)  !Implement floor to Ri
               RIQ=MIN(RIQ,RIQMAX)  !Implement ceiling to Ri
               
               RIQprint=RIQ !debug vertical profiles
               
               Knot = 
     &CP1*CP2*SURFVEL**2*EXP(-HP(L)*(1.0-ZZ(K))/CP1/SURFVEL) !K0(Z) = gamma*delta*omega^2*exp(-Z/gamma/omega)
               SFAB=Knot/(1.0+CP3*RIQ) !K(Z) = K0(Z)/[1+sigma*Ri(Z)]
             ELSEIF(HP(L)*(1.0-ZZ(K))<=ZMET(L))THEN !From Thermocline to Metalimnion: Vertical Eddy Diffusivity Computation  (Ztherm > Z >= Zmet)
!Calculation of K(Ztherm): KZTH
               Knottherm=
     &CP1*CP2*SURFVEL**2*EXP(-ZTHERM(L)/CP1/SURFVEL) !K0(Ztherm)
               BVFZtherm=
     &MAX(
     &GP/(1.0+B(L,KZTHERM(L)))*
     &ABS((B(L,KZTHERM(L)+1)-B(L,KZTHERM(L))))*HPI(L)*DZIG(KZTHERM(L)),
     &1.0E-18) !Brunt Vaisala frequency at Ztherm: BVF(Ztherm)
               RIQZtherm=
     &BVFZtherm*(CP1*CP2)**2*
     &EXP(2.*ZTHERM(L)/CP1/SURFVEL) !Ri(Ztherm)
               RIQZtherm=MAX(RIQZtherm,RIQMIN)  !Implement floor of Ri
               RIQZtherm=MIN(RIQZtherm,RIQMAX)  !Implement ceiling of Ri
               
               
               RIQprint=RIQZtherm !debug vertical profiles

               KZTH=Knottherm/(1.0+CP3*RIQZtherm) !K(Ztherm)
               IF(BVF>0.0)THEN
                 SFAB=alpha*KZTH*(BVFZtherm/BVF)**P2 !K(z) for Ztherm > Z >= Zmet
               ELSE
                 SFAB=0.0 !Avoid divide by zero errors
               ENDIF
             ELSE !In hypolimnion Zmet > Z > Zmax (or HP(L)
!Calculation of K(Zmet): KZM
               Knotmet=
     &CP1*CP2*SURFVEL**2*EXP(-ZMET(L)/CP1/SURFVEL)
!               BVFZmet=
!     &MAX(
!     &GP/(1.0+B(L,KZMET(L)))*
!     &ABS((B(L,KZMET(L)+1)-B(L,KZMET(L))))*HPI(L)*DZIG(KZMET(L)),
!     &    1.0E-6) !Brunt Vaisala frequency at Zmet: BFV(Zmet)
               BVFZmet=
     &MAX(-GP/(1.0E3*(1.0+B(L,KZMET(L))))*1E-5,
     &    1.0E-18) !Brunt Vaisala frequency at Zmet: BFV(Zmet) Using density gradient threshold 1E-5
               RIQZmet=
     &BVFZmet*(CP1*CP2)**2*EXP(2.*ZMET(L)/CP1/SURFVEL) !Ri(Zmet)
               RIQZmet=MAX(RIQZmet,RIQMIN)  !Implement floor of Ri
               RIQZmet=MIN(RIQZmet,RIQMAX)  !Implement ceiling of Ri
                              
               RIQprint=RIQZmet !debug vertical profiles
               

               KZM=Knotmet/(1.0+CP3*RIQZmet) !K(Zmet)
               IF(BVF>0.0)THEN
                 SFAB= 
     &(HP(L)-HP(L)*(1.0-ZZ(K)))/(HP(L)-ZMET(L))*
     &KZM*(BVFZmet/BVF)**P2 !K(z) for Zmet > Z >= Zmax (or HP(L)) 
               ELSE
                 SFAB=0.0
               ENDIF
             ENDIF
!Following Pacanowski and Philander (1981), used by Rao et al. (2004). 
!Need Richardson number, RIQ. Use the default calculation from EFDC (Brunt-Vaisala frequency divided by vertical shear).
!Vertical shear: [(del u)/(del z)]^2 + [(del v)/(del z)]^2
              SH2=((U(L,K+1)-U(L,K))*HPI(L)*DZIG(K))**2+
     &            ((V(L,K+1)-V(L,K))*HPI(L)*DZIG(K))**2
              RIQ=BVF/MIN(SH2+1.0E-6,RIQMAX)
!Richardson number is Brunt-Vaisala frequency divided by vertical shear 
!Previously calculated SFAV in EFDC following Mellor-Yamada (1982).
!SFAV=0.3920*(1.+8.6736*RIQ)/((1.+30.192*RIQ)*(1.+6.1272*RIQ))
!AVTMP=AVCON*SFAV*DML(L,K)*HP(L)*QQSQR(L,K)+AVO 
 
!New formulation by Pacanowski and Philander (1981)
!Kv0 is an adjustable parameter:  Kv0 =1.0 E-2 m2/s  
!AVO is background eddy diffusivity value (user input)
              SFAV = 0.01*(1.0+5.0*RIQ)**(-2)
!Include on/off switch (AVCON), add backgrounds, apply limits, and then harmonically smooth        
              ABTMP=AVCON*SFAB+ABO
              AVTMP=AVCON*SFAV+AVO
              AVMAX=MAX(AVMAX,AVTMP)
              ABMAX=MAX(ABMAX,ABTMP)
              AVMIN=MIN(AVMIN,AVTMP)
              ABMIN=MIN(ABMIN,ABTMP)
              AV(L,K)=       SQRT(AV(L,K)*AVTMP*HPI(L)) !Harmonic (time average) time smoothing
              AB(L,K)=SCB(L)*SQRT(AB(L,K)*ABTMP*HPI(L)) !Harmonic (time average) time smoothing
                            
            ENDDO
          ENDIF
        ENDDO  
      ENDIF
      ! *** NOW APPLY MAXIMUM, IF REQURIED
      IF(ISAVBMX.GE.1)THEN  
!$OMP PARALLEL PRIVATE(AVTMP,ABTMP)
        DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
          DO L=2,LA  
            AVTMP=AVMX*HPI(L)  
            ABTMP=ABMX*HPI(L)  
            AV(L,K)=MIN(AV(L,K),AVTMP)  
            AB(L,K)=MIN(AB(L,K),ABTMP)  
          ENDDO  
!$OMP END DO NOWAIT
        ENDDO  
!$OMP END PARALLEL
      ENDIF  
!$OMP PARALLEL PRIVATE(LS)
      DO K=1,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=2,LA  
          LS=LSC(L)  
c pmc          AVUI(L,K)=2./(AV(L,K)+AV(L-1,K))  
c pmc          AVVI(L,K)=2./(AV(L,K)+AV(LS,K))  
cferg          AVUI(L,K)=(1.+SUB(L))/(AV(L,K)+SUB(L)*AV(L-1,K))
cferg          AVVI(L,K)=(1.+SVB(L))/(AV(L,K)+SVB(L)*AV(LS,K))
               AVUI(L,K) = 1./AV(L,K)
               AVVI(L,K) = 1./AV(L,K)
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
      DO K=2,KS  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
        DO L=2,LA  
          AQ(L,K)=0.205*(AV(L,K-1)+AV(L,K))  
        ENDDO  
!$OMP END DO NOWAIT
      ENDDO  
!$OMP DO SCHEDULE(STATIC,CHUNKSIZE)
      DO L=2,LA  
        AQ(L,1)=0.205*AV(L,1)  
        AQ(L,KC)=0.205*AV(L,KS)  
      ENDDO  
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN  
      END  

