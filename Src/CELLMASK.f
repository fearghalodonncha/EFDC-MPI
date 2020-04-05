      SUBROUTINE CELLMASK  
C  
C CHANGE RECORD  
C **  SUBROUTINE CELLMASK CONVERTS LAND CELLS TO WATER CELLS BY  
C **  MASKING VARIABLES.  DEPTHS IN THE MASKED CELLS SHOULD BE INPUT AT  
C **  THE END OF THE DXDY.INP FILE.  
C  
      USE GLOBAL  
            LOGICAL    CELL_INSIDE_DOMAIN
      OPEN(1,FILE='MASK.INP',STATUS='UNKNOWN')  
      DO NS=1,6  
        READ(1,1)  
      ENDDO  
    1 FORMAT(80X)  
      READ(1,*,ERR=1000) MMASK  
      DO M=1,MMASK  
        READ(1,*,ERR=1000) II,JJ,MTYPE
        I = XLOC(II)
        J = YLOC(JJ)
        L=LIJ(I,J)  
        IF (CELL_INSIDE_DOMAIN(L)) THEN
            IF(MTYPE.EQ.1)THEN
              SUB(L)=0.
              SUBO(L)=0.
              UHDYE(L)=0.
              UHDY1E(L)=0.
              UHDY2E(L)=0.
              DO K=1,KC
                U(L,K)=0.
                U1(L,K)=0.
                U2(L,K)=0.
                UHDY(L,K)=0.
                UHDY1(L,K)=0.
                UHDY2(L,K)=0.
              ENDDO
            ENDIF
            IF(MTYPE.EQ.2)THEN
              SVB(L)=0.
              SVBO(L)=0.
              VHDXE(L)=0.
              VHDX1E(L)=0.
              VHDX2E(L)=0.
              DO K=1,KC
                V(L,K)=0.
                V1(L,K)=0.
                V2(L,K)=0.
                VHDX(L,K)=0.
                VHDX1(L,K)=0.
                VHDX2(L,K)=0.
              ENDDO
            ENDIF
            IF(MTYPE.EQ.3)THEN   ! *** PMC
              ! *** U Face
              SUB(L)=0.
              SUBO(L)=0.
              UHDYE(L)=0.
              UHDY1E(L)=0.
              UHDY2E(L)=0.
              DO K=1,KC
                U(L,K)=0.
                U1(L,K)=0.
                U2(L,K)=0.
                UHDY(L,K)=0.
                UHDY1(L,K)=0.
                UHDY2(L,K)=0.
              ENDDO
              ! *** V Face
              SVB(L)=0.
              SVBO(L)=0.
              VHDXE(L)=0.
              VHDX1E(L)=0.
              VHDX2E(L)=0.
              DO K=1,KC
                V(L,K)=0.
                V1(L,K)=0.
                V2(L,K)=0.
                VHDX(L,K)=0.
                VHDX1(L,K)=0.
                VHDX2(L,K)=0.
              ENDDO
            ENDIF
            IF(MTYPE.EQ.4)THEN   ! *** PMC - Change to MTYPE 4 for isolated waters
              LN=LNC(L)
              SUB(L)=0.
              SUBO(L)=0.
              UHDYE(L)=0.
              UHDY1E(L)=0.
              UHDY2E(L)=0.
              SUB(LEAST(L))=0.
              SUBO(LEAST(L))=0.
              UHDYE(LEAST(L))=0.
              UHDY1E(LEAST(L))=0.
              UHDY2E(LEAST(L))=0.
              SVB(L)=0.
              SVBO(L)=0.
              VHDXE(L)=0.
              VHDX1E(L)=0.
              VHDX2E(L)=0.
              SVB(LN)=0.
              SVBO(LN)=0.
              VHDXE(LN)=0.
              VHDX1E(LN)=0.
              VHDX2E(LN)=0.
              P(L)=0.
              P1(L)=0.
              DO K=1,KC
                B(L,K)=0.
                B1(L,K)=0.
                SAL(L,K)=0.
                SAL1(L,K)=0.
                TEM(L,K)=TEMO
                TEM1(L,K)=TEMO
                DYE(L,K)=0.
                DYE1(L,K)=0.
                SED(L,K,1)=0.
                SED1(L,K,1)=0.
                QQ(L,K)=0.
                QQ1(L,K)=0.
                QQL(L,K)=0.
                QQL1(L,K)=0.
                U(L,K)=0.
                U1(L,K)=0.
                U2(L,K)=0.
                UHDY(L,K)=0.
                UHDY1(L,K)=0.
                UHDY2(L,K)=0.
                U(LEAST(L),K)=0.
                U1(LEAST(L),K)=0.
                U2(LEAST(L),K)=0.
                UHDY(LEAST(L),K)=0.
                UHDY1(LEAST(L),K)=0.
                UHDY2(LEAST(L),K)=0.
                V(L,K)=0.
                V1(L,K)=0.
                V2(L,K)=0.
                VHDX(L,K)=0.
                VHDX1(L,K)=0.
                VHDX2(L,K)=0.
                V(LN,K)=0.
                V1(LN,K)=0.
                V2(LN,K)=0.
                VHDX(LN,K)=0.
                VHDX1(LN,K)=0.
                VHDX2(LN,K)=0.
              ENDDO
            ENDIF
      END IF
      ENDDO  
      CLOSE(1)  
C  
C **  WRITE READ ERRORS ON CELLMASK  
C  
      GOTO 1002  
 1000 WRITE(6,1001)  
 1001 FORMAT('  READ ERROR ON FILE MASK.INP ')  
      STOP  
 1002 CONTINUE  
      RETURN  
      END  

