      SUBROUTINE SURFPLT  
C  
C CHANGE RECORD  
C **  SUBROUTINE SURFPLT WRITES FILES TO CONTOUR FREE SURFACE  
C **  ELEVATION  
C  
      USE GLOBAL
!      CHARACTER*80 TITLE  
      INTEGER*8 TIME
      LOGICAL,SAVE:: FIRSTTIME=.FALSE.
C  
C *** EE BEGIN BLOCK  
C  
      INTEGER*4    VER  
C  
C *** EE END BLOCK  
C  
      IF(IPPHXY.LE.2)THEN  
        IF(JSPPH.NE.1) GOTO 300  
        OPEN(10,FILE='SURFCON'//ans(partid2)//'.OUT')  
        CLOSE(10,STATUS='DELETE')  
        JSPPH=0
        FIRSTTIME = .TRUE.
  300   CONTINUE  
        IF(ISDYNSTP.EQ.0)THEN  
          TIME=DT*FLOAT(N)+TCON*TBEGIN  
        ELSE  
          TIME=TIMESEC  
        ENDIF  
        OPEN(10,FILE='SURFCON'//ans(partid2)//'.OUT',POSITION='APPEND')  
        WRITE (10,100)N,TIME,PARTID,LA,XPAR(1),XPAR(IC),YPAR(1),YPAR(JC)  
        IF(IPPHXY.EQ.0)THEN  
          DO L=2,LA  
            SURFEL=BELV(L)+HP(L)  
            WRITE(10,201)SURFEL,BELV(L),HP(L),  
     &            HBED(L,KBT(L)),HBEDA(L)  
          ENDDO  
        ENDIF  
        IF(IPPHXY.EQ.1)THEN  
          DO L=2,LA  
            I = XPAR(IL(L))
            J = YPAR(JL(L))
            SURFEL=BELV(L)+HP(L)  
            WRITE(10,200)I,J,SURFEL,BELV(L),HP(L),  
     &            HBED(L,KBT(L)),HBEDA(L)  
          ENDDO  
        ENDIF  
        IF(IPPHXY.EQ.2)THEN  
          DO L=2,LA  
            SURFEL=BELV(L)+HP(L)  
            WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),SURFEL,BELV(L)
     &                 ,HP(L),HBED(L,KBT(L)),HBEDA(L)  
          ENDDO  
        ENDIF  
        CLOSE(10)  
      ENDIF  
C  
C *** EE BEGIN BLOCK  
C *** OUTPUT EFDC EXPLORER FORMAT.  DO NOT CHANGE OUTPUTS!  
C ***                               MUST EXACTLY MATCH EFDC_EXPLORER INP  
C  
      IF(IPPHXY.EQ.3)THEN  
        LINES=LA-1  
        IF(JSPPH.EQ.1)THEN  
          OPEN(10,FILE='EE_WS.OUT',STATUS='UNKNOWN')
          CLOSE(10,STATUS='DELETE')  
          OPEN(10,FILE='EE_WS.OUT',STATUS='UNKNOWN',  
     &        ACCESS='SEQUENTIAL',FORM='UNFORMATTED')  
          VER=101  
          WRITE(10)VER,IC,JC,LINES  
          JSPPH=0  
          CLOSE(10)  
        ENDIF  
        OPEN(10,FILE='EE_WS.OUT',POSITION='APPEND',STATUS='OLD',  
     &      FORM='UNFORMATTED')  !UNFORMATTED
        IF(ISDYNSTP.EQ.0)THEN  
          TIME=DT*FLOAT(N)+TCON*TBEGIN  
        ELSE  
          TIME=TIMESEC  
        ENDIF  
        TIME=TIME/86400.  
        IF(ISDYNSTP.EQ.0)THEN  
          DELT=DT  
        ELSE  
          DELT=DTDYN  
        ENDIF  
        WRITE (10)N,TIME,DELT  
        DO L=2,LA  
          WRITE(10)HP(L)  
        ENDDO 
 !       CALL FLUSH(10) !This was commented out by Fearghal
        CLOSE(10)  
      ENDIF  
C  
C *** EE END BLOCK  
C  
   99 FORMAT(A80)  
  100 FORMAT(2I14,I3,5I7)  
  101 FORMAT(2I10)  
  200 FORMAT(2I5,1X,9E14.5)  
  201 FORMAT(9E14.5)  
  250 FORMAT(12E12.4)  
      RETURN  
      END  

