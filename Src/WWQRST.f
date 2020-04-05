      SUBROUTINE WWQRST  
C  
C CHANGE RECORD  
C WRITE SPATIAL DISTRIBUTIONS AT THE END OF SIMULATION TO UNIT IWQORST.  
C  
      USE GLOBAL  
C  
C WRITE ASCII RESTART FILE:  
C  
      OPEN(1,FILE='WQWCRST.OUT',STATUS='UNKNOWN')  
      CLOSE(1,STATUS='DELETE')  
      OPEN(1,FILE='WQWCRST.OUT',STATUS='UNKNOWN')  
      IF(ISDYNSTP.EQ.0)THEN  
        TIME=DT*FLOAT(N)+TCON*TBEGIN  
        TIME=TIME/TCON  
      ELSE  
        TIME=TIMESEC/TCON  
      ENDIF  
      WRITE(1,101) N,TIME  
      WRITE(1,102)  
      NWQV0=NWQV  
      IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1  
      DO L=2,LA  
        DO K=1,KC  
          WRITE(1,90) L,K,(WQV(L,K,NW),NW=1,NWQV0)  
        ENDDO  
      ENDDO  
      CLOSE(1)  
C  
C ALSO WRITE BINARY RESTART FILE:  
C  
   90 FORMAT(2I5, 1P, 23E12.4)  
  101 FORMAT('CC  WQ RESTART FILE TIME STEP, TIME = ',I10,F12.5)  
  102 FORMAT('C   L    K  BC          BD          BG          ',  
     &    'RPOC        LPOC        DOC         ',  
     &    'RPOP        LPOP        DOP         PTO4        ',  
     &    'RPON        LPON        DON         AMN         ',  
     &    'NIT         SU          SA          COD         ',  
     &    'DO          TAM         FCB        MALG')  
      RETURN  
      END  

