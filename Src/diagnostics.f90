SUBROUTINE DGN_PUV  
        USE GLOBAL
        OPEN(1,FILE='FUV.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
        WRITE(1,1001)N,ISTL  
        DO L=2,LA  
          WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HU(L),P1(L), &
             P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)  
        ENDDO  
        CLOSE(1)  
        IF(N.EQ.1)THEN  
          OPEN(1,FILE='FUV1.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,1001)N,ISTL  
          DO L=2,LA  
            WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HU(L),P1(L), &
               P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
        IF(N.EQ.2)THEN  
          OPEN(1,FILE='FUV2.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,1001)N,ISTL  
          DO L=2,LA  
            WRITE(1,1001)IL(L),JL(L),UHDY1E(L),HRUO(L),HU(L),P1(L), &
              P1(L-1),TSX1(L),TBX1(L),FCAXE(L),FPGXE(L),FXE(L)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
 1001 FORMAT(2I5,10(1X,E12.4))  
END


SUBROUTINE DGN_PUV_2
        USE GLOBAL
        OPEN(1,FILE='FP.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
        WRITE(1,1001)N,ISTL  
        DO L=2,LA  
          WRITE(1,1001)IL(L),JL(L),FP(L),FUHDYE(L),FUHDYE(L+1), &  
          FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)  
        ENDDO  
        CLOSE(1)  
        IF(N.EQ.1)THEN  
          OPEN(1,FILE='FP.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,1001)N,ISTL  
          DO L=2,LA  
            WRITE(1,1001)IL(L),JL(L),FP(L),FUHDYE(L),FUHDYE(L+1), &
            FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
        IF(N.EQ.2)THEN  
          OPEN(1,FILE='FP2.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,1001)N,ISTL  
          DO L=2,LA  
            WRITE(1,1001)IL(L),JL(L),FP(L),FUHDYE(L),FUHDYE(L+1), &  
            FVHDXE(L),FVHDXE(LNC(L)),QSUME(L),RIFTR(L),EVAPSW(L)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
 1001 FORMAT(2I5,10(1X,E12.4))  
END


SUBROUTINE DGN_PUV_3
        USE GLOBAL
        OPEN(1,FILE='EQCOEF.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
        WRITE(1,1001)N,ISTL  
        DO L=2,LA  
          SURFTMP=GI*P(L)  
          WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L), & 
           FP(L),SURFTMP  
        ENDDO  
        CLOSE(1)  
        IF(N.EQ.1)THEN  
          OPEN(1,FILE='EQCOEF1.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,1001)N,ISTL  
          DO L=2,LA  
            SURFTMP=GI*P(L)  
            WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L), &  
             FP(L),SURFTMP  
          ENDDO  
          CLOSE(1)  
        ENDIF  
        IF(N.EQ.2)THEN  
          OPEN(1,FILE='EQCOEF2.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,1001)N,ISTL  
          DO L=2,LA  
            SURFTMP=GI*P(L)  
            WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L), &  
             FP(L),SURFTMP  
          ENDDO  
          CLOSE(1)  
        ENDIF  

        OPEN(1,FILE='EQTERM.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
        WRITE(1,1001)N,ISTL  
        DO L=2,LA  
          WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L), & 
           HRVO(L),HU(L),HV(L)  
        ENDDO  
        CLOSE(1)  
        IF(N.EQ.1)THEN  
          OPEN(1,FILE='EQTERM1.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,1001)N,ISTL  
          DO L=2,LA  
            WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L), & 
            HRVO(L),HU(L),HV(L)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
        IF(N.EQ.2)THEN  
          OPEN(1,FILE='EQTERM2.OUT',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,1001)N,ISTL  
          DO L=2,LA  
            WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L), & 
               HRVO(L),HU(L),HV(L)  
          ENDDO  
          CLOSE(1)  
        ENDIF  
 1001 FORMAT(2I5,10(1X,E12.4))  
END
