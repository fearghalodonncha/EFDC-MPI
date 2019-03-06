
program readdata_testing

implicit none

integer :: ULOC, LA_ZCAL,LA_PRAN,LA_HORDIF,LA_BEGTI, LA_ENDTI, LA_FREQ
integer :: NPD

OPEN(ULOC,FILE='DRIFTER.INP',ACTION='READ')
  CALL READSTR(ULOC)
  READ(ULOC,*) LA_ZCAL,LA_PRAN,LA_HORDIF !05 MAY 2009: NEW STRUCTURE 
  
	write(*,*) LA_ZCAL,LA_PRAN,LA_HORDIF
	write(*,*) 'testing..........'
  
  CALL READSTR(ULOC)
  READ(ULOC,*) LA_BEGTI, LA_ENDTI, LA_FREQ            !UPDATED 23-04-09
  !LA_FREQ = LA_FREQ/1440.                             !Output Frequency 
  
	print *, LA_BEGTI, LA_ENDTI, LA_FREQ
	write(*,*) 'testing..........'
  
  CALL READSTR(ULOC)
  READ(ULOC,*) NPD

	print *, NPD
	write(*,*) 'testing..........'
end program readdata_testing

SUBROUTINE READSTR(UINP)   !******************************************************************
  INTEGER(4),INTENT(IN)::UINP
  CHARACTER(200)::STR
  DO WHILE (1==1)
    READ(UINP,'(A)') STR
    STR=ADJUSTL(STR)
    IF (STR(1:1).NE.'*') THEN
      BACKSPACE(UINP)
      RETURN
    ENDIF
  ENDDO
END SUBROUTINE
