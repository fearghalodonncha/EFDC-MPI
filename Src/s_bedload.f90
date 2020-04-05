SUBROUTINE BEDLOADJ
  USE GLOBAL     
  IMPLICIT NONE 
  INTEGER::I,J,L,K,LW
  !PT: real values are written in DOUBLE PRECISION 7/16/08.
  DOUBLE PRECISION,DIMENSION(LCM)::VELMAG
  DOUBLE PRECISION::UTMP,VTMP 
  !
  !  University of California, Santa Barbara
  !  Craig Jones and Wilbert Lick
  !
  !  Bedload transport subroutine based on Van Rijn's transport
  !  Equations.
  !
  !  REVISION DATE :  May 24, 2006
  !  Craig Jones and Scott James
  !***************************************************************
  !  Calculate Percentage of erosion into suspension PSUS
  !  and whether the cell has bedload or not BLFLAG
  FORALL(L=2:LA)
     FORALL(K=1:NSCM)
        USW(L,K)=SQRT(TAU(L))/DWS(K) !DWS is settling speed.  USW is the shear velocity.
     ENDFORALL
  ENDFORALL
  !Loop below determines if bedload exists or not.  There are three regimes of tranport in this loop.
  !the first conditional in the where check for a large enough particle's diameter and small enough
  !shear velocity.  If the particle is too small or the shear velocity is too large, then all the sediment 
  !transport is in the suspended load, specified by suspended probability (PSUS = 1).  If the particle
  !is large enough and the shear velocity is small enough then we have two situations. In the first case, 
  !shear stress tau is smaller than the critical shear velocity or if shear velocity is negative or zero
  !then there is neither bedload transport or suspended load transport.  Otherwise, both bedload and suspended
  !load transport exists.  Also calculated is the probability of suspension for suspended load PSUS (eqn. 8).  
  LCM_LOOP_1:DO L=2,LA
     WHERE(D50(1:NSCM)>=200.0.AND.USW(L,1:NSCM)<4.0)
        WHERE(TAU(L)<=TCRE(1:NSCM).OR.USW(L,1:NSCM)<=0.0)
           BLFLAG(L,1:NSCM)=0
           PSUS(L,1:NSCM)=0.0
        ELSEWHERE                              
           BLFLAG(L,1:NSCM)=1
           PSUS(L,1:NSCM)=MAX((LOG(USW(L,1:NSCM))-LOG(SQRT(TCRSUS(1:NSCM))/ &
           DWS(1:NSCM)))/(LOG(4.0)-LOG(SQRT(TCRSUS(1:NSCM))/DWS(1:NSCM))),0.d0)
        ENDWHERE
     ELSEWHERE
        BLFLAG(L,1:NSCM)=0
        PSUS(L,1:NSCM)=1.0
     ENDWHERE
  ENDDO LCM_LOOP_1
  SED_LOOP:DO K=1,NSCM
     LCM_LOOP_2:DO L=2,LA
        VELMAG(L)=SQRT(U(L,1)**2+V(L,1)**2)
        TRANS(L,K)=MAX((TAU(L)-TCRE(K))/TCRE(K),0.d0) !eqn. 21
        DZBL(L,K)=MIN(100.d0*HP(L),D50(K)/10000.d0*0.3d0*DISTAR(K)**0.7d0*SQRT(TRANS(L,K))) !(eqn. 20b) don't allow bedload height to exceed water column depth
        IF(DZBL(L,K)/=0.0.AND.DZBL_LAST(L,K)/=0.0)CBL(1,L,K)=CBL(1,L,K)*DZBL_LAST(L,K)/DZBL(L,K)
        BLVEL(L,K)=1.5*TRANS(L,K)**0.6*SQRT(((SEDDENS/WATERDENS) -1.0)*980.0*D50(K)/10000.0)      !eqn. 20a
     ENDDO LCM_LOOP_2
     WHERE(VELMAG(2:LA)>0.0)
        UBL(2:LA,K)=BLVEL(2:LA,K)*U(2:LA,1)/VELMAG(2:LA)
        VBL(2:LA,K)=BLVEL(2:LA,K)*V(2:LA,1)/VELMAG(2:LA)
     ELSEWHERE           
        UBL(2:LA,K)=0.0
        VBL(2:LA,K)=0.0
     ENDWHERE
     IF(ISSLOPE /=0)THEN !if bedslope is calculated
       DO L=2,LA
         UTMP=UBL(L,K) !save original x-bedload velocity
         VTMP=VBL(L,K) !save original x-bedload velocity
         UBL(L,K)=ALPHA_PX(L)*UBL(L,K) !modify by pitch angle
         VBL(L,K)=ALPHA_PY(L)*VBL(L,K) !modify by roll angle
         IF(UBL(L,K)>VBL(L,K))THEN !find dominant velocity direction
           VBL(L,K)=VBL(L,K)+ALPHA_RX(L,K)*UTMP !Bedload velocity (x is dominant roll rirection)
           UBL(L,K)=UBL(L,K)-ALPHA_RX(L,K)*VBL(L,K) !as impacted by
           UBL(L,K)=UBL(L,K)+ALPHA_RY(L,K)*VTMP !bedslope; see (secondary roll due to y)
           VBL(L,K)=VBL(L,K)-ALPHA_RY(L,K)*UBL(L,K) !Lesser (2004) Ikeda (1982)
         ELSE
           UBL(L,K)=UBL(L,K)+ALPHA_RY(L,K)*VTMP !Bedload velocity (y is dominant roll direction)
           VBL(L,K)=VBL(L,K)-ALPHA_RY(L,K)*UBL(L,K) !as impacted by
           VBL(L,K)=VBL(L,K)+ALPHA_RX(L,K)*UTMP !bedslope; see (secondary roll due to x)
           UBL(L,K)=UBL(L,K)-ALPHA_RX(L,K)*VBL(L,K) !Lesser (2004) Ikeda (1982)
         ENDIF
       ENDDO
     ENDIF
     
     !**********************************************************************!
     ! All the equations below are solving the pde in eqn.18.
     ! X Bedload flux at I-1/2 interface
     LCM_LOOP_3:DO L=2,LA
        IF(IL(LWEST(L))==0)THEN
          LW=0
        ELSE
          LW=LIJ(IL(LWEST(L)),JL(L))
        ENDIF
        IF(LW==0)THEN !is it a western boundary?
           XBLFLUX(L,K)=DT/(DXU(L)*100.0)*CBL(1,L,K)*UBL(L,K)*DZBL(L,K) !no directionality
        ELSE !it is an internal cell (or even on the east boundary)
           IF(UBL(L,K)>=0.0)THEN !check which direction is upstream
              XBLFLUX(L,K)=DT/(DXU(L)*100.0)*CBL(1,LW,K)*UBL(L,K)*0.5*(DZBL(LW,K)+DZBL(L,K)) !use bedload concentration from the east because flow is east to west
           ELSE
              XBLFLUX(L,K)=DT/(DXU(L)*100.0)*CBL(1,L ,K)*UBL(L,K)*0.5*(DZBL(LW,K)+DZBL(L,K)) !use bedload concentration from the currect cell to calculate flux across the east face of the current cell
           ENDIF
        ENDIF
     ENDDO LCM_LOOP_3
     
     ! Y Bedload flux at J-1/2 interface
     LCM_LOOP_5:DO L=2,LA
        IF(LSC(L)==0)THEN !is it a southern boundary?
          YBLFLUX(L,K)=DT/(DYV(L)*100.0)*CBL(1,L,K)*DZBL(L,K) !no directionality
        ELSE !it is an internal cell (or even on the north boundary)
          IF(VBL(L,K)>=0)THEN !check which direction is upstream
            YBLFLUX(L,K)=DT/(DYV(L)*100.0)*CBL(1,LSC(L),K)*VBL(L,K)*0.5*(DZBL(LSC(L),K)+DZBL(L,K)) !use bedload concentration from the south because flow is south to north
          ELSE
            YBLFLUX(L,K)=DT/(DYV(L)*100.0)*CBL(1,L     ,K)*VBL(L,K)*0.5*(DZBL(LSC(L),K)+DZBL(L,K)) !use bedload concentration from the current cell to calculate flux across the south face of the current cell
          ENDIF
        ENDIF
     ENDDO LCM_LOOP_5
     
     !**********************************************************************!
     !  Transport Equation Interior Elements
     ! 
     LCM_LOOP_4:DO L=2,LA
        I=IL(L)
        J=JL(L)
        IF(DZBL(L,K)>0.0)THEN
          IF(LIJ(I+1,J)>0.AND.LIJ(I,J+1)>0)THEN !is it an interior cell?
            CBL(2,L,K)=CBL(1,L,K)+(XBLFLUX(L,K)-XBLFLUX(LIJ(I+1,J),K)+YBLFLUX(L,K)-YBLFLUX(LIJ(I,J+1),K)+QBSED(L,K))/DZBL(L,K) !full bedload flux calculation
          ELSEIF(LIJ(I+1,J)>0.AND.LIJ(I,J+1)==0)THEN !is it a northern boundary?
            CBL(2,L,K)=CBL(1,L,K)+(XBLFLUX(L,K)-XBLFLUX(LIJ(I+1,J),K)+QBSED(L,K))/DZBL(L,K) !no y bedload flux component
          ELSEIF(LIJ(I+1,J)==0.AND.LIJ(I,J+1)>0)THEN !is it an eastern boundary?
            CBL(2,L,K)=CBL(1,L,K)+(YBLFLUX(L,K)-YBLFLUX(LIJ(I,J+1),K)+QBSED(L,K))/DZBL(L,K) !no x bedload flux component
          ELSE !is it the northeast corner?
            CBL(2,L,K)=CBL(1,L,K)+QBSED(L,K)/DZBL(L,K) !no x or y bedload flux components
          ENDIF
        ENDIF
        IF(DZBL(L,K)==0.OR.CBL(2,L,K)<0)CBL(2,L,K)=0.0 
     ENDDO LCM_LOOP_4
  ENDDO SED_LOOP
  CBL(1,2:LA,1:NSCM)=CBL(2,2:LA,1:NSCM)
  
  !CQBEDLOADX(:,:)=XBLFLUX(:,:)
  !CQBEDLOADY(:,:)=YBLFLUX(:,:)
  DZBL_LAST(:,:)=DZBL(:,:)
  RETURN 
END SUBROUTINE BEDLOADJ
