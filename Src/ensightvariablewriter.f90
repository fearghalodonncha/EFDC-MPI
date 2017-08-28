SUBROUTINE ENSIGHTVARIABLEWRITE(ENSNAME,ENSVARIABLE,IUNIT,LENTEMP)
USE GLOBAL
IMPLICIT NONE
INTEGER, INTENT(IN):: IUNIT, LENTEMP
REAL, INTENT(IN), DIMENSION(LENTEMP) :: ENSVARIABLE
INTEGER::L
CHARACTER(20)::ENSNAME					!VARIABLE TAG

 !**********************************************************
 !Writing the variable file.
  WRITE(IUNIT,*)'BEGIN TIME STEP'
 !First we write the description of the variable (variable name).
  WRITE(IUNIT,*)TRIM(ENSNAME)
 !part here should coincide with the part written in .geo file.  We'll let it equal 1 always since our grid doesn't change.
  WRITE(IUNIT,*)'part'
  WRITE(IUNIT,11) 1
  WRITE(IUNIT,*)'block'			!ensight formatting requirements
  11 FORMAT(I3)			
  12 FORMAT(3I3)
  13 FORMAT(E12.5)
  !data must be output in the same format as the .geo file. 
  DO L=2,LENTEMP
        WRITE(IUNIT,13) ENSVARIABLE(L)
  ENDDO
  WRITE(IUNIT,*)'END TIME STEP'
RETURN
END SUBROUTINE ENSIGHTVARIABLEWRITE