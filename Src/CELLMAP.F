      SUBROUTINE CELLMAP  
C  
C **  SUBROUTINE CELLMAP GENERATES CELL MAPPINGS  
C CHANGE RECORD  
C  
      USE GLOBAL  
      INTEGER,ALLOCATABLE,DIMENSION(:)::IL_GL
      INTEGER,ALLOCATABLE,DIMENSION(:)::JL_GL
      ALLOCATE(IL_GL(LC_GLOBAL+2))
      ALLOCATE(JL_GL(LC_GLOBAL+2))
#ifdef key_mpi
             CALL CHILDGRID
#else
      DO I = 1,IC
         XPAR(I) = I
      END DO
      DO J = 1,JC
         YPAR(J) = J
      END DO
#endif
      DO J =1,JC
        DO I =1,IC   
          IJCT(I,J)=IJCT_GLOBAL(XPAR(I),YPAR(J))
          IJCTLT(I,J)=IJCTLT_GLOBAL(XPAR(I),YPAR(J))
        END DO
      END DO
      IF (PARTID .EQ. 0) THEN
        OPEN(37,FILE='CELLMAP1.DAT',status='unknown')
        L=2
        DO J =1,JC_GLOBAL
          DO I =1,IC_GLOBAL
            IF(IJCT_GLOBAL(i,j).EQ.5)THEN
              IF(PARTID==1) write(37,*)I,J,L,JC_GLOBAL,IC_GLOBAL,LC_GLOBAL
              L = L + 1
            END IF
          ENDDO
        ENDDO
      END IF
      LIJ(:,:) = 0
C **  SET 1D CELL INDEX SEQUENCE AND MAPPINGS  
      L=2
      DO J =1,JC
        DO I =1,IC
          IF(IJCT(I,J).EQ.5)THEN         
            IL_GL(L)=XPAR(I)
            JL_GL(L)=YPAR(J)
            IL(L) = I
            JL(L) = J
            LCT(L)=IJCT_GLOBAL(IL_GL(L),JL_GL(L))
            LIJ(I,J)=L
            IF(PARTID==MASTER_TASK) write(37,*)I,J,IL_GL(L),JL_GL(L),LIJ(I,J),L,JC,IC
            L = L + 1
          END IF
        ENDDO
      ENDDO
      close(37)
      LA=L-1  
      LCTT=L  
      IF(LCTT.NE.LC)THEN  
        LC = LCTT
      ENDIF  
      IL(1)=0  
      IL(LC)=0  
      JL(1)=0  
      JL(LC)=0  
c      WRITE(1,601)LA  
      WRITE(7,601)LA  
c      CLOSE(1)  
  601 FORMAT('  LA=',I10,//)  
 1616 FORMAT(2I10)  
 1617 FORMAT(' LC =',I5,' DETERMINED IN CELLMAP',  
     &    ' INCONSISTENT WITH INPUT VALUE =',I5//)  

      DO I =1, IC
        DO J =1,JC
          IF (XPAR(I).LE.0) IJCT(I,J) = 0
          IF (XPAR(I).GT.IC_GLOBAL) IJCT(I,J) = 0
          IF (YPAR(J).LE.0) IJCT(I,J) = 0
          IF (YPAR(J).GT.JC_GLOBAL)  IJCT(I,J) = 0
        END DO
      END DO

!! extract coordinates of computational footprint only (i.e. excluding ghost zones) 
!! for in particular conjugate gradient computation
      LL = 0
      DO I = 3,IC-2
        DO J = 3,JC-2
          IF (IJCT(I,J) == 5) THEN
            LL = LL + 1
          END IF
        END DO
      END DO
      ALLOCATE(L_CONG(LL))
      CONGDOM = LL
      LL = 0
      L_CONG = 0
      DO I = 3,IC-2
        DO J = 3,JC-2
          IF (IJCT(I,J) == 5) THEN
            LL = LL + 1
            L_CONG(LL) = LIJ(I,J)
          END IF
        END DO
      END DO


      open (3,FILE='CELLK'//ANS(PARTID2)//'.DAT',STATUS='UNKNOWN')
      write(3,*)'IGRID=',XPAR(1),XPAR(IC),PARTID
      DO J =JC,1,-1
        write(3,33)YPAR(J),(IJCT(I,J),I=1,IC)
      END DO
      CLOSE(3)
33    FORMAT(I4,2X,3000I1)

      IF(ISWASP.LE.5)THEN  
        L=2  
        DO J=1,JC  
          DO I=1,IC  
            IF(IJCTLT(I,J).GT.0.AND.IJCTLT(I,J).LT.9)THEN  
              ILLT(L)=I  
              JLLT(L)=J  
              LCTLT(L)=IJCTLT(I,J)  
              LIJLT(I,J)=L  
              L=L+1  
            ENDIF  
          ENDDO  
        ENDDO  
        LALT=L-1  
        LCLT=L  
      ENDIF  
      IF(ISWASP.GE.6)THEN  
        L=2  
        DO J=1,JC  
          DO I=1,IC  
            IF(IJCTLT(I,J).GT.0.AND.IJCTLT(I,J).LT.7)THEN  
              ILLT(L)=I  
              JLLT(L)=J  
              LCTLT(L)=IJCTLT(I,J)  
              LIJLT(I,J)=L  
              L=L+1  
            ENDIF  
          ENDDO  
        ENDDO  
        LALT=L-1  
        LCLT=L  
      ENDIF  
      WRITE(7,1616)LALT,LCLT
  
C **  ASSIGN RED AND BLACK CELL SEQUENCES  (PMC - NOT FUNCTIONAL)
C  
      !IF(IRVEC.NE.9)THEN  
      !ENDIF  
C  
C **  SET NORTH AND SOUTH CELL IDENTIFIER ARRAYS  
C  
      LNC(1)=LC  
      LSC(1)=LC  
      LNEC(1)=LC  
      LNWC(1)=LC  
      LSEC(1)=LC  
      LSWC(1)=LC  
      LEAST(1)=LC
      LWEST(1)=LC
      LNC(LC)=1  
      LSC(LC)=1  
      LNEC(LC)=1  
      LNWC(LC)=1  
      LSEC(LC)=1  
      LSWC(LC)=1 
      LEAST(LC)=1
      LWEST(LC)=1 
      DO L=2,LA  
        I=IL(L)  
        J=JL(L)  
        IF(IJCT(I,J+1).EQ.9)THEN  
          LNC(L)=LC  
        ELSE  
          LNC(L)=LIJ(I,J+1)  
        ENDIF  
        IF(IJCT(I+1,J+1).EQ.9)THEN  
          LNEC(L)=LC  
        ELSE  
          LNEC(L)=LIJ(I+1,J+1)  
        ENDIF  
        IF(IJCT(I-1,J+1).EQ.9)THEN  
          LNWC(L)=LC  
        ELSE  
          LNWC(L)=LIJ(I-1,J+1)  
        ENDIF  
        IF(IJCT(I-1,J).EQ.9)THEN
          LWEST(L)=LC
        ELSE
          LWEST(L)=MAX(LIJ(I-1,J),1)
        ENDIF
        IF(IJCT(I+1,J).EQ.9)THEN
          LEAST(L)=LC
        ELSE
          LEAST(L)=MAX(LIJ(I+1,J),1)
        ENDIF
        IF(J-1<1)CYCLE
        IF(IJCT(I+1,J-1).EQ.9)THEN  
          LSEC(L)=LC  
        ELSE  
          LSEC(L)=LIJ(I+1,J-1)  
        ENDIF 
        IF(IJCT(I-1,J-1).EQ.9)THEN  
          LSWC(L)=LC  
        ELSE  
          LSWC(L)=LIJ(I-1,J-1)  
        ENDIF  
        IF(IJCT(I,J-1).EQ.9)THEN  
          LSC(L)=LC  
        ELSE  
          LSC(L)=LIJ(I,J-1)  
        ENDIF  
      ENDDO
C  
C **  MODIFY NORTH-SOUTH CELL MAPPING FOR PERIOD GRID IN N-S DIRECTION  
C  
      IF(ISPGNS.GE.1)THEN  
        DO NPN=1,NPNSBP  
          LS=LIJ(ISPNS(NPN),JSPNS(NPN))  
          LSC(LS)=LIJ(INPNS(NPN),JNPNS(NPN))  
          IF( IJCT(INPNS(NPN)+1,JNPNS(NPN)).EQ.9)THEN  
            LSEC(LS)=LC  
          ELSE  
            LSEC(LS)=LIJ(INPNS(NPN)+1,JNPNS(NPN))  
          ENDIF  
          IF( IJCT(INPNS(NPN)-1,JNPNS(NPN)).EQ.9)THEN  
            LSWC(LS)=LC  
          ELSE  
            LSWC(LS)=LIJ(INPNS(NPN)-1,JNPNS(NPN))  
          ENDIF  
          LN=LIJ(INPNS(NPN),JNPNS(NPN))  
          LNC(LN)=LIJ(ISPNS(NPN),JSPNS(NPN))  
          IF( IJCT(ISPNS(NPN)+1,JSPNS(NPN)).EQ.9)THEN  
            LNEC(LN)=LC  
          ELSE  
            LNEC(LN)=LIJ(ISPNS(NPN)+1,JSPNS(NPN))  
          ENDIF  
          IF( IJCT(ISPNS(NPN)-1,JSPNS(NPN)).EQ.9)THEN  
            LNWC(LN)=LC  
          ELSE  
            LNWC(LN)=LIJ(ISPNS(NPN)-1,JSPNS(NPN))  
          ENDIF  
        ENDDO  
      ENDIF  
C  
C **  SET LT NORTH AND SOUTH CELL IDENTIFIER ARRAYS  
C  
      LNCLT(1)=LCLT  
      LSCLT(1)=LCLT  
      LNCLT(LC)=1  
      LSCLT(LC)=1  
      DO L=2,LALT  
        I=ILLT(L)  
        J=JLLT(L)
        IF(IJCTLT(I,J+1).EQ.9)THEN  
          LNCLT(L)=LCLT  
        ELSE  
          LNCLT(L)=LIJLT(I,J+1)  
        ENDIF  
        IF(J-1<1)CYCLE
        IF(IJCTLT(I,J-1).EQ.9)THEN  
          LSCLT(L)=LCLT  
        ELSE  
          LSCLT(L)=LIJLT(I,J-1)  
        ENDIF  
      ENDDO

C
C  hnr start modification may 2002---------------------------------C
C
C **  MODIFY lt NORTH-SOUTH CELL MAPPING FOR PERIOD GRID IN N-S DIRECTION
C
      IF (ISPGNS.GE.1) THEN
        DO NPN=1,NPNSBP
C
C     SET NORTH CELL SOUTH OF SOUTH CELL
C
          LS=LIJLT(ISPNS(NPN),JSPNS(NPN))
          LSCLT(LS)=LIJLT(INPNS(NPN),JNPNS(NPN))
C
C     SET SOUTH CELL NORTH OF NORTH CELL
C
          LN=LIJLT(INPNS(NPN),JNPNS(NPN))
          LNCLT(LN)=LIJLT(ISPNS(NPN),JSPNS(NPN))

        END DO
      END IF
C
C  hnr end modification may 2002 -----------------------------------C
C  
C **  SET NORTH, SOUTH, EAST AND WEST RED CELL IDENTIFIER ARRAYS  (NOT FUNCTIONAL)
C  
      !IF(IRVEC.NE.9)THEN  
      !ENDIF  
C  
C **  SET NORTH, SOUTH, EAST AND WEST BLACK CELL IDENTIFIER ARRAYS  (NOT FUNCTIONAL)
C  
      !IF(IRVEC.NE.9)THEN  
      !ENDIF  
C  
C **  DIAGNOSE OF RED-BLACK CELL MAPPING  
C  
      !IF(IRVEC.EQ.2.OR.IRVEC.EQ.9) GOTO 220  
C  ! PMC - NOT USED
C **  RED CELL LOOP  (PMC - RED/BLACK ITERATIVE SOLVER NOT FUNCIONAL)
C  
      !DO LR=1,NRC  
      !ENDDO  
C  
C **  BLACK CELL LOOP  
C  
      !DO LB=1,NBC  
      !ENDDO  
      
!  220 CONTINUE  
  101 FORMAT(' LR,LTMP = ',2I6/)  
  102 FORMAT(' LR,LTMP = ',2I6/)  
  103 FORMAT(' LN= 1, LR,LTMP,ITMP,JTMP = ',4I6/)  
  104 FORMAT(' LN=LC, LR,LTMP,ITMP,JTMP = ',4I6/)  
  105 FORMAT(' NERR, LR,LTMP,ITMP,JTMP,INTMP,JNTMP = ',6I6/)  
  106 FORMAT(' NERR, LR,LTMP,ITMP,JTMP,INTMP,JNTMP = ',6I6/)  
  107 FORMAT(' LS= 1, LR,LTMP,ITMP,JTMP = ',4I6/)  
  108 FORMAT(' LS=LC, LR,LTMP,ITMP,JTMP = ',4I6/)  
  109 FORMAT(' SERR, LR,LTMP,ITMP,JTMP,ISTMP,JSTMP = ',6I6/)  
  110 FORMAT(' SERR, LR,LTMP,ITMP,JTMP,ISTMP,JSTMP = ',6I6/)  
  111 FORMAT(' LE= 1, LR,LTMP,ITMP,JTMP = ',4I6/)  
  112 FORMAT(' LE=LC, LR,LTMP,ITMP,JTMP = ',4I6/)  
  113 FORMAT(' EERR, LR,LTMP,ITMP,JTMP,IETMP,JETMP = ',6I6/)  
  114 FORMAT(' EERR, LR,LTMP,ITMP,JTMP,IETMP,JETMP = ',6I6/)  
  115 FORMAT(' LW= 1, LR,LTMP,ITMP,JTMP = ',4I6/)  
  116 FORMAT(' LW=LC, LR,LTMP,ITMP,JTMP = ',4I6/)  
  117 FORMAT(' WERR, LR,LTMP,ITMP,JTMP,IWTMP,JWTMP = ',6I6/)  
  118 FORMAT(' WERR, LR,LTMP,ITMP,JTMP,IWTMP,JWTMP = ',6I6/)  
  201 FORMAT(' LB,LTMP = ',2I6/)  
  202 FORMAT(' LB,LTMP = ',2I6/)  
  203 FORMAT(' LN= 1, LB,LTMP,ITMP,JTMP = ',4I6/)  
  204 FORMAT(' LN=LC, LB,LTMP,ITMP,JTMP = ',4I6/)  
  205 FORMAT(' NERR, LB,LTMP,ITMP,JTMP,INTMP,JNTMP = ',6I6/)  
  206 FORMAT(' NERR, LB,LTMP,ITMP,JTMP,INTMP,JNTMP = ',6I6/)  
  207 FORMAT(' LS= 1, LB,LTMP,ITMP,JTMP = ',4I6/)  
  208 FORMAT(' LS=LC, LB,LTMP,ITMP,JTMP = ',4I6/)  
  209 FORMAT(' SERR, LB,LTMP,ITMP,JTMP,ISTMP,JSTMP = ',6I6/)  
  210 FORMAT(' SERR, LB,LTMP,ITMP,JTMP,ISTMP,JSTMP = ',6I6/)  
  211 FORMAT(' LE= 1, LB,LTMP,ITMP,JTMP = ',4I6/)  
  212 FORMAT(' LE=LC, LB,LTMP,ITMP,JTMP = ',4I6/)  
  213 FORMAT(' EERR, LB,LTMP,ITMP,JTMP,IETMP,JETMP = ',6I6/)  
  214 FORMAT(' EERR, LB,LTMP,ITMP,JTMP,IETMP,JETMP = ',6I6/)  
  215 FORMAT(' LW= 1, LB,LTMP,ITMP,JTMP = ',4I6/)  
  216 FORMAT(' LW=LC, LB,LTMP,ITMP,JTMP = ',4I6/)  
  217 FORMAT(' WERR, LB,LTMP,ITMP,JTMP,IWTMP,JWTMP = ',6I6/)  
  218 FORMAT(' WERR, LB,LTMP,ITMP,JTMP,IWTMP,JWTMP = ',6I6/)  
  119 FORMAT('       SUB(LETMP) = ',F10.2/)  
  120 FORMAT('       SUB(LTMP)  = ',F10.2/)  
C  
C **  DEFINE MAPPING TO 3D GRAPHICS GRID  
C  
      IF(ISCLO.EQ.0.OR.NWGG.EQ.0)THEN  
        IG=IC  
        JG=JC  
      ELSE  
        OPEN(1,FILE='GCELLMP.INP',STATUS='UNKNOWN')  
        READ(1,1111)  
        READ(1,1111)  
        READ(1,1111)  
        READ(1,1111)  
        READ(1,*)IG,JG  
        DO NW=1,NWGG  
          READ(1,*)IGTMP,JGTMP,ICOMP,JCOMP  
          LTMP=LIJ(ICOMP,JCOMP)  
          IWGG(NW)=IGTMP  
          JWGG(NW)=JGTMP  
          LWGG(NW)=LTMP  
        ENDDO  
        CLOSE(1)  
      ENDIF  
 1111 FORMAT(80X)  
      DEALLOCATE(IL_GL)
      DEALLOCATE(JL_GL)

      RETURN  
      END  

