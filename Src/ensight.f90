SUBROUTINE ENSIGHT
  USE GLOBAL
  IMPLICIT NONE
  !This file writes the .geo and .U (where U is the variable of interest)
  !for ensight.  The .geo and .U files will later be read in by the
  !file .case for the ensight format.
  INTEGER:: I, J, K, LENTEMP
  REAL:: RHOw, HTCAP, HTCONT, HTEMP
  CHARACTER(20)::ENSCHARTMP
  REAL,DIMENSION(KC)::CTEMP1
  REAL, ALLOCATABLE, DIMENSION(:) :: ENSTEMP1, ENSTEMP2, ENSTEMP3, ENSTEMP4, ENSTEMP5
  REAL, ALLOCATABLE, DIMENSION(:) :: ENSTEMP6, ENSTEMP7, ENSTEMP8, ENSTEMP9, ENSTEMP10
  REAL, ALLOCATABLE, DIMENSION(:) :: ENSTEMP11, ENSTEMP12, ENSTEMP13, ENSTEMP14, ENSTEMP15
  REAL, ALLOCATABLE, DIMENSION(:) :: ENSTEMP16, ENSTEMP17, ENSTEMP18, ENSTEMP19, ENSTEMP20
  REAL, ALLOCATABLE, DIMENSION(:) :: ENSTEMP21, ENSTEMP22, ENSTEMP23, ENSTEMP24, AVGSED, UTMPS, VTMPS
  LENTEMP=MAXVAL(LIJ(3:IC-2,3:JC-2))
  ALLOCATE(ENSTEMP1(LENTEMP))
  ALLOCATE(ENSTEMP2(LENTEMP))
  ALLOCATE(ENSTEMP3(LENTEMP))
  ALLOCATE(ENSTEMP4(LENTEMP))
  ALLOCATE(ENSTEMP5(LENTEMP))
  ALLOCATE(ENSTEMP6(LENTEMP))
  ALLOCATE(ENSTEMP7(LENTEMP))
  ALLOCATE(ENSTEMP8(LENTEMP))
  ALLOCATE(ENSTEMP9(LENTEMP))
  ALLOCATE(ENSTEMP10(LENTEMP))
  ALLOCATE(ENSTEMP11(LENTEMP))
  ALLOCATE(ENSTEMP12(LENTEMP))
  ALLOCATE(ENSTEMP13(LENTEMP))
  ALLOCATE(ENSTEMP14(LENTEMP))
  ALLOCATE(ENSTEMP15(LENTEMP))
  ALLOCATE(ENSTEMP16(LENTEMP))
  ALLOCATE(ENSTEMP17(LENTEMP))
  ALLOCATE(ENSTEMP18(LENTEMP))
  ALLOCATE(ENSTEMP19(LENTEMP))
  ALLOCATE(ENSTEMP20(LENTEMP))
  ALLOCATE(ENSTEMP21(LENTEMP))
  ALLOCATE(ENSTEMP22(LENTEMP))
  ALLOCATE(ENSTEMP23(LENTEMP))
  ALLOCATE(ENSTEMP24(LENTEMP))
  ALLOCATE(AVGSED(LENTEMP))
  ALLOCATE(UTMPS(LENTEMP))
  ALLOCATE(VTMPS(LENTEMP))
    !**********************************************************
  !Writing the geometry file for each time step.
  WRITE(114,*)'BEGIN TIME STEP'
  WRITE(114,*)'EnSight Model Geometry File'
  WRITE(114,*)'EnSight 8.2.5'
  WRITE(114,*)'node id given'
  WRITE(114,*)'element id assign'
  WRITE(114,*)'extents'
  WRITE(114,10) MINVAL(DLON(2:MAXVAL(LIJ(3:IC-2,3:JC-2)))),MAXVAL(DLON(2:MAXVAL(LIJ(3:IC-2,3:JC-2))))
  WRITE(114,10) MINVAL(DLAT(2:MAXVAL(LIJ(3:IC-2,3:JC-2)))),MAXVAL(DLAT(2:MAXVAL(LIJ(3:IC-2,3:JC-2))))
  WRITE(114,10) 0.0, 0.0
  WRITE(114,*) 'part'
  WRITE(114,11) 1
  WRITE(114,11) 1
  WRITE(114,*)'block'
  WRITE(114,12) IC-4, JC-4, 1
  DO J=3,JC-2
    DO I=3,IC-2
        WRITE(114,13) DLON(LIJ(I,J))
    ENDDO
  ENDDO
  DO J=3,JC-2
    DO I=3,IC-2
        WRITE(114,13) DLAT(LIJ(I,J))
    ENDDO
  ENDDO
  DO J=3,JC-2
    DO I=3,IC-2
        WRITE(114,13) 0.0
    ENDDO
  ENDDO
  WRITE(114,*)'node_ids'
  DO J=3,JC-2
    DO I=3,IC-2
        WRITE(114,14) LIJ(I,J)-1
    ENDDO
  ENDDO
  10 FORMAT(2E12.5)
  11 FORMAT(I3)				
  12 FORMAT(3(6X,I3))		
  13 FORMAT(E13.6)			!VB changed frm e12.5
  14 FORMAT(6X,I4) !ensight is 'picky' so we have to put in the 6X!!!!
  WRITE(114,*)'END TIME STEP'
  !**************************************************************
!PT writing out ensight.case header.
      CLOSE(113,STATUS='DELETE')
      OPEN(UNIT=113,FILE='ensight.case',FORM='FORMATTED')
      VARUNCOUNT=VARUNCOUNT+1
      WRITE(113,*)'#'
      WRITE(113,*)'#'
      WRITE(113,*)'#'
      WRITE(113,*)'# Case File: ensight.case' 
      WRITE(113,*)'FORMAT'  
      WRITE(113,*)'type: ensight gold'
      WRITE(113,*)'GEOMETRY'
      WRITE(113,*)'model:    1           1       ensight.geo'
      WRITE(113,*)
      WRITE(113,*)'VARIABLE'
      IF (ENSIGHT1 .GT. 0) THEN
        WRITE(113,*)'scalar per node: 1     1      U ensight.U' 
      ENDIF
      IF (ENSIGHT2 .GT. 0) THEN
        WRITE(113,*)'scalar per node: 1     1      V ensight.V'
      ENDIF
sedloop: IF(ISTRAN(6) .GT. 0) THEN
            IF (ENSIGHT3 .GT. 0) THEN
                WRITE(113,*)'scalar per node: 1     1      TAU ensight.TAU' 
            ENDIF
            IF (ENSIGHT4 .GT. 0) THEN
                WRITE(113,*)'scalar per node: 1     1      D50 ensight.D50' 
            ENDIF
            IF (ENSIGHT5 .GT. 0) THEN
                WRITE(113,*)'scalar per node: 1     1      CBL ensight.CBL' 
            ENDIF
            IF (ENSIGHT6 .GT. 0) THEN
                WRITE(113,*)'scalar per node: 1     1      SED ensight.SED' 
            ENDIF
      ENDIF sedloop
wqloop:      IF (ISTRAN(8) .GT. 0) THEN
        IF (ENSIGHT7 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      CHC ensight.CHC' 
        ENDIF
        IF (ENSIGHT8 .GT. 0) THEN
			 WRITE(113,*)'scalar per node: 1     1      CHD ensight.CHD' 	
        ENDIF  
        IF (ENSIGHT9 .GT. 0) THEN
			 WRITE(113,*)'scalar per node: 1     1      CHG ensight.CHG' 
        ENDIF   
        IF (ENSIGHT10 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      ROC ensight.ROC' 
        ENDIF  
        IF (ENSIGHT11 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      LOC ensight.LOC' 
        ENDIF
        IF (ENSIGHT12 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      DOC ensight.DOC' 
        ENDIF
        IF (ENSIGHT13 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      ROP ensight.ROP' 
        ENDIF
        IF (ENSIGHT14 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      LOP ensight.LOP' 
        ENDIF
        IF (ENSIGHT15 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      DOP ensight.DOP' 
        ENDIF
        IF (ENSIGHT16 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      P4D ensight.P4D' 
        ENDIF
        IF (ENSIGHT17 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      RON ensight.RON' 
        ENDIF
        IF (ENSIGHT18 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      LON ensight.LON' 
        ENDIF
        IF (ENSIGHT19 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      DON ensight.DON' 
        ENDIF
        IF (ENSIGHT20 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      NHX ensight.NHX' 
        ENDIF
        IF (ENSIGHT21 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      NOX ensight.NOX' 
        ENDIF
        IF (ENSIGHT22 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      SUU ensight.SUU' 
        ENDIF
        IF (ENSIGHT23 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      SAA ensight.SAA' 
        ENDIF
        IF (ENSIGHT24 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      COD ensight.COD' 
        ENDIF
        IF (ENSIGHT25 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      DOX ensight.DOX' 
        ENDIF
        IF (ENSIGHT26 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      TAM ensight.TAM' 
        ENDIF
        IF (ENSIGHT27 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      FCB ensight.FCB' 
        ENDIF
        IF (ENSIGHT28 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      CO2 ensight.CO2' 
        ENDIF
        IF (ENSIGHT29 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      MAC ensight.MAC' 
        ENDIF
        IF (ENSIGHT30 .GT. 0) THEN
            WRITE(113,*)'scalar per node: 1     1      HEAT ensight.HEAT' 
        ENDIF
      ENDIF wqloop
      IF ((ENSIGHT31 .GT. 0) .AND.(ISTRAN(2) .GT. 0)) THEN
            WRITE(113,*)'scalar per node: 1     1      TEMP ensight.TEMP'
      ENDIF
      WRITE(113,*)
      WRITE(113,*)'TIME'
      WRITE(113,*)'time set: 1'
      WRITE(113,*)'number of steps: ', VARUNCOUNT
      WRITE(113,*)'time values:'
      DO K=1,VARUNCOUNT
        WRITE(113,13) TBEGIN+FLOAT(K-1)*DT*FLOAT(ISHPRT-1)/86400		!VB CONVERSION TO DAYS
      ENDDO
      WRITE(113,*)'FILE'
      WRITE(113,*)'file set: 1'
      WRITE(113,*)'number of steps: ', VARUNCOUNT
      !The rest of the ensight.case file will be written in the ensight.f90 subroutine since there are 
      !time step information that is only known later.
!************************************************************ 
!***** variable calculations ***********
    !VB HEAT CONTENT VARIABLE	
	RHOw=1000.	   !KG/M3
	HTCAP=4187.	   !J/KG-K
	HTCONT=0.0						!VB INITIALIZED
	DO I=3,IC-2
	    DO J=3,JC-2
	        HTEMP=RHOw*HTCAP*SUM((TEM(LIJ(I,J),1:KC)+273.15)*DZC(1:KC))*DXYP(LIJ(I,J))
	        HTCONT=HTCONT+HTEMP   !HEAT CONTENT IN JOULES
	    ENDDO
	ENDDO
    DO I=3,IC-2
        DO J=3,JC-2
            UTMPS(LIJ(I,J))     = SUM((U(LIJ(I,J),1:KC)*CUE(LIJ(I,J))+V(LIJ(I,J),1:KC)*CVE(LIJ(I,J)))*DZC(1:KC))
            VTMPS(LIJ(I,J))     = SUM((U(LIJ(I,J),1:KC)*CUN(LIJ(I,J))+V(LIJ(I,J),1:KC)*CVN(LIJ(I,J)))*DZC(1:KC))
        IF (ISTRAN(6) .GT. 0) THEN
            CBLTOT(LIJ(I,J))    = 10.0*SUM(CBL(1,LIJ(I,J),1:NSED)*DZBL(LIJ(I,J),1:NSED))*DXYP(LIJ(I,J)) !g/cm^3*cm*m^2*(0.001*100*100)
            CTEMP1(1:KC)        = 0.001*SUM(SED(LIJ(I,J),1:KC,1:NSED))
            AVGSED(LIJ(I,J))    = SUM(CTEMP1(1:KC)*DZC(1:KC))*DXYP(LIJ(I,J))*HP(LIJ(I,J))
        ENDIF					
        IF (ISTRAN(8) .GT. 0) THEN
                ENSTEMP1(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,1))/FLOAT(KC)
                ENSTEMP2(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,2))/FLOAT(KC)
                ENSTEMP3(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,3))/FLOAT(KC)
                ENSTEMP4(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,4))/FLOAT(KC)
                ENSTEMP5(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,5))/FLOAT(KC)
                ENSTEMP6(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,6))/FLOAT(KC)
                ENSTEMP7(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,7))/FLOAT(KC)
                ENSTEMP8(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,8))/FLOAT(KC)
                ENSTEMP9(LIJ(I,J))  = SUM(WQV(LIJ(I,J),1:KC,9))/FLOAT(KC)
                ENSTEMP10(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,10))/FLOAT(KC)
                ENSTEMP11(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,11))/FLOAT(KC)
                ENSTEMP12(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,12))/FLOAT(KC)
                ENSTEMP13(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,13))/FLOAT(KC)
                ENSTEMP14(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,14))/FLOAT(KC)
                ENSTEMP15(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,15))/FLOAT(KC)
                ENSTEMP16(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,16))/FLOAT(KC)
                ENSTEMP17(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,17))/FLOAT(KC)
                ENSTEMP18(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,18))/FLOAT(KC)
                ENSTEMP19(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,19))/FLOAT(KC)
                ENSTEMP20(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,20))/FLOAT(KC)
                ENSTEMP21(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,21))/FLOAT(KC)
                ENSTEMP22(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,22))/FLOAT(KC)
                ENSTEMP23(LIJ(I,J)) = SUM(WQV(LIJ(I,J),1:KC,23))/FLOAT(KC)
                ENSTEMP24(LIJ(I,J)) = HTCONT
         ENDIF
         IF (ISTRAN(2) .GT. 0) THEN
            TWATER(LIJ(I,J))=SUM(TEM(LIJ(I,J),1:KC)*DZC(1:KC))
         ENDIF
        ENDDO
    ENDDO

!***** end of variable calculations ****   
! Now let's input the variable files.
    IF (ENSIGHT1 .GT. 0) THEN
        ENSCHARTMP = "U"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,UTMPS,118,LENTEMP)
    ENDIF
    IF (ENSIGHT2 .GT. 0) THEN
        ENSCHARTMP = "V"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,VTMPS,119,LENTEMP)
    ENDIF
IF (ISTRAN(6) .GT. 0 ) THEN
    IF (ENSIGHT3 .GT. 0) THEN
        ENSCHARTMP = "TAU"
        !CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,TAU,120,LENTEMP)
    ENDIF
    IF (ENSIGHT4 .GT. 0) THEN
        ENSCHARTMP = "D50"
        !CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,D50AVG,121,LENTEMP)
    ENDIF
    IF (ENSIGHT5 .GT. 0) THEN
        ENSCHARTMP = "CBL"
        !CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,CBLTOT,122,LENTEMP)
    ENDIF
    IF (ENSIGHT6 .GT. 0) THEN
        ENSCHARTMP = "SED"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,AVGSED,123,LENTEMP)
    ENDIF
ENDIF
IF (ISTRAN(8) .GT. 0) THEN
    IF (ENSIGHT7 .GT. 0) THEN
        ENSCHARTMP = "CHC"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP1,124,LENTEMP)
    ENDIF
    IF (ENSIGHT8 .GT. 0) THEN
		ENSCHARTMP = "CHD"				
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP2,125,LENTEMP)
    ENDIF
    IF (ENSIGHT9 .GT. 0) THEN
		ENSCHARTMP = "CHG"				
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP3,126,LENTEMP)
    ENDIF
    IF (ENSIGHT10 .GT. 0) THEN
        ENSCHARTMP = "ROC"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP4,127,LENTEMP)
    ENDIF
    IF (ENSIGHT11 .GT. 0) THEN
        ENSCHARTMP = "LOC"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP5,128,LENTEMP)
    ENDIF
    IF (ENSIGHT12 .GT. 0) THEN
        ENSCHARTMP = "DOC"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP6,129,LENTEMP)
    ENDIF
    IF (ENSIGHT13 .GT. 0) THEN
        ENSCHARTMP = "ROP"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP7,130,LENTEMP)
    ENDIF
    IF (ENSIGHT14 .GT. 0) THEN
        ENSCHARTMP = "LOP"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP8,131,LENTEMP)
    ENDIF
    IF (ENSIGHT15 .GT. 0) THEN
        ENSCHARTMP = "DON"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP9,132,LENTEMP)
    ENDIF
    IF (ENSIGHT16 .GT. 0) THEN
        ENSCHARTMP = "P4D"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP10,133,LENTEMP)
    ENDIF
    IF (ENSIGHT17 .GT. 0) THEN
        ENSCHARTMP = "RON"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP11,134,LENTEMP)
    ENDIF
    IF (ENSIGHT18 .GT. 0) THEN
        ENSCHARTMP = "LON"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP12,135,LENTEMP)
    ENDIF
    IF (ENSIGHT19 .GT. 0) THEN
        ENSCHARTMP = "DON"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP13,136,LENTEMP)
    ENDIF
    IF (ENSIGHT20 .GT. 0) THEN
        ENSCHARTMP = "NHX"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP14,137,LENTEMP)
    ENDIF
    IF (ENSIGHT21 .GT. 0) THEN
        ENSCHARTMP = "NOX"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP15,138,LENTEMP)
    ENDIF
    IF (ENSIGHT22 .GT. 0) THEN
        ENSCHARTMP = "SUU"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP16,139,LENTEMP)
    ENDIF
    IF (ENSIGHT23 .GT. 0) THEN
        ENSCHARTMP = "SAA"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP17,140,LENTEMP)
    ENDIF
    IF (ENSIGHT24 .GT. 0) THEN
        ENSCHARTMP = "COD"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP18,141,LENTEMP)
    ENDIF
    IF (ENSIGHT25 .GT. 0) THEN
        ENSCHARTMP = "DOX"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP19,142,LENTEMP)
    ENDIF
    IF (ENSIGHT26 .GT. 0) THEN
        ENSCHARTMP = "TAM"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP20,143,LENTEMP)
    ENDIF
    IF (ENSIGHT27 .GT. 0) THEN
        ENSCHARTMP = "FCB"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP21,144,LENTEMP)
    ENDIF
    IF (ENSIGHT28 .GT. 0) THEN
        ENSCHARTMP = "CO2"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP22,145,LENTEMP)
    ENDIF
    IF (ENSIGHT29 .GT. 0) THEN
        ENSCHARTMP = "MAC"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP23,146,LENTEMP)
    ENDIF
    IF (ENSIGHT30 .GT. 0) THEN
        ENSCHARTMP = "HEAT"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,ENSTEMP24,147,LENTEMP)
    ENDIF
ENDIF
IF (ISTRAN(2) .GT. 0) THEN
    IF(ENSIGHT31 .GT. 0) THEN
        ENSCHARTMP = "TEMP"
        CALL ENSIGHTVARIABLEWRITE(ENSCHARTMP,TWATER,148,LENTEMP)
    ENDIF
ENDIF
!************************************************************  
  	RETURN
END SUBROUTINE ENSIGHT