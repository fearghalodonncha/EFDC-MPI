      SUBROUTINE CALBAL2  
C  
C CHANGE RECORD  
C **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM,  
C **  AND ENERGY BALANCES  
C  
      USE GLOBAL  
	IMPLICIT NONE
	INTEGER::LL,K,LS,L,LN
C  
C **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES  
C  
      DO K=1,KC  
        DO LL=1,NCBS  
          L=LCBS(LL)  
          LN=LNC(L)  
          VOLOUT=VOLOUT-VHDX2(LN,K)*DZC(K)  
          SALOUT=SALOUT-MIN(VHDX2(LN,K),0.)*SAL1(LN,K)*DZC(K)  
     &        -MAX(VHDX2(LN,K),0.)*SAL1(L,K)*DZC(K)  
          DYEOUT=DYEOUT-MIN(VHDX2(LN,K),0.)*DYE1(LN,K)*DZC(K)  
     &        -MAX(VHDX2(LN,K),0.)*DYE1(L,K)*DZC(K)  
          PPEOUT=PPEOUT-VHDX2(LN,K)*G*DZC(K)*( 0.5*(BELV(L)+BELV(LN))  
     &        +0.125*(HP(L)+H2P(L)+HP(LN)+H2P(LN))*(Z(K)+Z(K-1)) )  
          BBEOUT=BBEOUT-MIN(VHDX2(LN,K),0.)*DZC(K)*GP*( BELV(LN)  
     &        +0.5*HP(LN)*(Z(K)+Z(K-1)) )*B1(LN,K)  
     &        -MAX(VHDX2(LN,K),0.)*DZC(K)*GP*( BELV(L)  
     &        +0.5*HP(L)*(Z(K)+Z(K-1)) )*B1(L,K)  
        ENDDO  
      ENDDO  
      DO K=1,KC  
        DO LL=1,NCBW  
          L=LCBW(LL)  
          VOLOUT=VOLOUT-UHDY2(LEAST(L),K)*DZC(K)  
          SALOUT=SALOUT-MIN(UHDY2(LEAST(L),K),0.)*SAL1(LEAST(L),K)*DZC(K)  
     &        -MAX(UHDY2(LEAST(L),K),0.)*SAL1(L,K)*DZC(K)  
          DYEOUT=DYEOUT-MIN(UHDY2(LEAST(L),K),0.)*DYE1(LEAST(L),K)*DZC(K)  
     &        -MAX(UHDY2(LEAST(L),K),0.)*DYE1(L,K)*DZC(K)  
          PPEOUT=PPEOUT-UHDY2(LEAST(L),K)*G*DZC(K)*( 0.5*(BELV(L)+BELV(LEAST(L)))  
     &        +0.125*(HP(L)+H2P(L)+HP(LEAST(L))+H2P(LEAST(L)))*(Z(K)+Z(K-1)) )  
          BBEOUT=BBEOUT-MIN(UHDY2(LEAST(L),K),0.)*DZC(K)*GP*( BELV(LEAST(L))  
     &        +0.5*HP(LEAST(L))*(Z(K)+Z(K-1)) )*B1(LEAST(L),K)  
     &        -MAX(UHDY2(LEAST(L),K),0.)*DZC(K)*GP*( BELV(L)  
     &        +0.5*HP(L)*(Z(K)+Z(K-1)) )*B1(L,K)  
        ENDDO  
      ENDDO  
      DO K=1,KC  
        DO LL=1,NCBE  
          L=LCBE(LL)  
          VOLOUT=VOLOUT+UHDY2(L,K)*DZC(K)  
          SALOUT=SALOUT+MIN(UHDY2(L,K),0.)*SAL1(L,K)*DZC(K)  
     &        +MAX(UHDY2(L,K),0.)*SAL1(LWEST(L),K)*DZC(K)  
          DYEOUT=DYEOUT+MIN(UHDY2(L,K),0.)*DYE1(L,K)*DZC(K)  
     &        +MAX(UHDY2(L,K),0.)*DYE1(LWEST(L),K)*DZC(K)  
          PPEOUT=PPEOUT+UHDY2(L,K)*G*DZC(K)*( 0.5*(BELV(L)+BELV(LWEST(L)))  
     &        +0.125*(HP(L)+H2P(L)+HP(LWEST(L))+H2P(LWEST(L)))*(Z(K)+Z(K-1)) )  
          BBEOUT=BBEOUT+MIN(UHDY2(L,K),0.)*DZC(K)*GP*(BELV(L)  
     &        +0.5*HP(L)*(Z(K)+Z(K-1)) )*B1(L,K)  
     &        +MAX(UHDY2(L,K),0.)*DZC(K)*GP*(BELV(LWEST(L))  
     &        +0.5*HP(LWEST(L))*(Z(K)+Z(K-1)) )*B1(LWEST(L),K)  
        ENDDO  
      ENDDO  
      DO K=1,KC  
        DO LL=1,NCBN  
          L=LCBN(LL)  
          LS=LSC(L)  
          VOLOUT=VOLOUT+VHDX2(L,K)*DZC(K)  
          SALOUT=SALOUT+MIN(VHDX2(L,K),0.)*SAL1(L,K)*DZC(K)  
     &        +MAX(VHDX2(L,K),0.)*SAL1(LS,K)*DZC(K)  
          DYEOUT=DYEOUT+MIN(VHDX2(L,K),0.)*DYE1(L,K)*DZC(K)  
     &        +MAX(VHDX2(L,K),0.)*DYE1(LS,K)*DZC(K)  
          PPEOUT=PPEOUT+VHDX2(L,K)*G*DZC(K)*( 0.5*(BELV(L)+BELV(LS))  
     &        +0.125*(HP(L)+H2P(L)+HP(LS)+H2P(LS))*(Z(K)+Z(K-1)) )  
          BBEOUT=BBEOUT+MIN(VHDX2(L,K),0.)*DZC(K)*GP*( BELV(L)  
     &        +0.5*HP(L)*(Z(K)+Z(K-1)) )*B1(L,K)  
     &        +MAX(VHDX2(L,K),0.)*DZC(K)*GP*( BELV(LS)  
     &        +0.5*HP(LS)*(Z(K)+Z(K-1)) )*B1(LS,K)  
        ENDDO  
      ENDDO  
      RETURN  
      END  

