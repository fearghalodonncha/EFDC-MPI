      SUBROUTINE TOXCHEM  
C  
C CHANGE RECORD  
C **  SUBROUTINE CALSND CALCULATES NONCOHESIVE SEDIMENT SETTLING,  
C **  DEPOSITION AND RESUSPENSION AND IS CALLED FOR SSEDTOX  
C  
      USE GLOBAL  
      IF(ISTRAN(5).GE.1)THEN  
        DO NT=1,NTOX  
C  
C **     NOTES:  
C        BULK DECAY COEFFICIENT  
C        VOLITIZATION  
C          VOLTOX(NT)  
C          RMOLTX(NT)=MOLECULAR WEIGHT  
C        PHOTOLOSIS  
C          RKTOXP(NT)=BASE RATE  
C          SKTOXP(NT)=SOLAR RADIATION AT BASE RATE  
C  
          DO K=1,KC  
            DO L=2,LA  
              CDECAYW(L,K)=1./(1.+DELT*RKTOXW(NT))  
            ENDDO  
          ENDDO  
          DO K=1,KC  
            DO L=2,LA  
              TOX(L,K,NT)=CDECAYW(L,K)*TOX(L,K,NT)  
            ENDDO  
          ENDDO  
          DO K=1,KB  
            DO L=2,LA  
              CDECAYB(L,K)=1./(1.+DELT*RKTOXB(NT))  
            ENDDO  
          ENDDO  
          DO K=1,KB  
            DO L=2,LA  
              TOXB(L,K,NT)=CDECAYB(L,K)*TOXB(L,K,NT)  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  
      RETURN  
      END  

