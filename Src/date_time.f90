
 FUNCTION IDAY(YYYY, MM, DD) RESULT(IVAL)
!====IDAY IS A COMPANION TO CALEND; GIVEN A CALENDAR DATE, YYYY, MM,
!           DD, IDAY IS RETURNED AS THE DAY OF THE YEAR.
!           EXAMPLE: IDAY(1984, 4, 22) = 113

 INTEGER, INTENT(IN) :: YYYY, MM, DD
 INTEGER             :: IVAL

 IVAL = 3055*(MM+2)/100 - (MM+10)/13*2 -91 +  &
       (1-(MODULO(YYYY, 4)+3)/4              &
        + (MODULO(YYYY, 100) + 99)/100 -        &
       (MODULO(YYYY, 400)+399)/400)*(MM+10)/13 + DD

 RETURN
 END FUNCTION IDAY


FUNCTION JD_OUT(YYYY, MM, DD) RESULT(IVAL)
  INTEGER, INTENT(IN)  :: YYYY
  INTEGER, INTENT(IN)  :: MM
  INTEGER, INTENT(IN)  :: DD
  INTEGER*8            :: IVAL
!              DATE ROUTINE JD(YYYY, MM, DD) CONVERTS CALENDER DATE TO
!              JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
!              EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
!    EXAMPLE JD(1970, 1, 1) = 2440588
  IVAL = DD - 32075 + 1461*(YYYY+4800+(MM-14)/12)/4 +  &
        367*(MM-2-((MM-14)/12)*12)/12 - 3*((YYYY+4900+(MM-14)/12)/100)/4
  RETURN
 END FUNCTION JD_OUT

 SUBROUTINE CDATE(JD,YYYY,MM,DD)
!====GIVEN A JULIAN DAY NUMBER, NNNNNNNN, YYYY,MM,DD ARE RETURNED AS THE
!         CALENDAR DATE. JD = NNNNNNNN IS THE JULIAN DATE FROM AN EPOCH
!         IN THE VERY DISTANT PAST.  SEE CACM 1968 11(10):657,
!         LETTER TO THE EDITOR BY FLIEGEL AND VAN FLANDERN.
!    EXAMPLE CALL CDATE(2440588, YYYY, MM, DD) RETURNS 1970 1 1 .

 INTEGER,INTENT(IN)  ::  JD
 INTEGER,INTENT(OUT) ::  YYYY,MM,DD
 INTEGER L

 L = JD + 68569
 N = 4*L/146097
 L = L - (146097*N + 3)/4
 YYYY = 4000*(L+1)/1461001
 L = L - 1461*YYYY/4 + 31
 MM = 80*L/2447
 DD = L - 2447*MM/80
 L = MM/11
 MM = MM + 2 - 12*L
 YYYY = 100*(N-49) + YYYY + L
 RETURN
 END SUBROUTINE CDATE




      FUNCTION REF_ATT (R_TIME, R_DATE)
!
!=======================================================================
!                                                                      !
!  This function encodes the relative time attribute that gives the    !
!  elapsed interval since a specified reference time.  The "units"     !
!  attribute takes the form "time-unit since reference-time".          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     r_time     Time-reference (real; %Y%m%d.%f, for example,         !
!                  20020115.5 for 15 Jan 2002, 12:0:0).                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ref_att    Time-reference for "units" attribute (string).        !
!     r_date     Calendar date vector (real):                          !
!                  r_date(1) => reference date (yyyymmdd.f).           !
!                  r_date(2) => year.                                  !
!                  r_date(3) => year day.                              !
!                  r_date(4) => month.                                 !
!                  r_date(5) => day.                                   !
!                  r_date(6) => hour.                                  !
!                  r_date(7) => minute.                                !
!                  r_date(8) => second.                                !
!                                                                      !
!=======================================================================
!
!
      
!
!  Imported variable declarations.
!
      USE GLOBAL
      REAL(WP), INTENT(IN) :: R_TIME
      REAL(WP), DIMENSION(8), INTENT(OUT) :: R_DATE
      CHARACTER (LEN=19) :: REF_ATT
!
!  Local variable declarations.
!
      INTEGER(WP) :: IDAY, IHOUR, ISEC, IYEAR, LEAP, MINUTE, MONTH
      INTEGER(WP), DIMENSION(13) :: IYD =                                   &
     &         (/ 1,32,60,91,121,152,182,213,244,274,305,335,366 /)
      INTEGER(WP), DIMENSION(13) :: IYDL =                                  &
     &         (/ 1,32,61,92,122,153,183,214,245,275,306,336,367 /)
      REAL(WP) :: DAY, SEC, YDAY
      CHARACTER (LEN=19) :: TEXT
!
!-----------------------------------------------------------------------
!  Decode reference time.
!-----------------------------------------------------------------------
!
      IYEAR=MAX(1,INT(R_TIME*0.0001))
      MONTH=MIN(12,MAX(1,INT((R_TIME-REAL(IYEAR*10000))*0.01)))
      DAY=R_TIME-AINT(R_TIME*0.01)*100.0
      IDAY=INT(DAY)
      SEC=(DAY-AINT(DAY))*86400.0
      IHOUR=INT(SEC/3600.0)
      MINUTE=INT(MOD(SEC,3600.0)/60.0)
      ISEC=INT(MOD(SEC,60.0))
!
!-----------------------------------------------------------------------
!  Get year day.
!-----------------------------------------------------------------------
!
      LEAP=MOD(IYEAR,4)
      IF (LEAP.EQ.0) THEN
        YDAY=REAL(IYDL(MONTH))+REAL(IDAY)-1.0
      ELSE
        YDAY=REAL(IYD(MONTH))+REAL(IDAY)-1.0
      END IF
!
!-----------------------------------------------------------------------
!  Build output date vector.
!-----------------------------------------------------------------------
!
      R_DATE(1)=R_TIME
      R_DATE(2)=REAL(IYEAR)
      R_DATE(3)=MAX(1.0,YDAY)
      R_DATE(4)=REAL(MONTH)
      R_DATE(5)=MAX(1.0,REAL(IDAY))
      R_DATE(6)=REAL(IHOUR)
      R_DATE(7)=REAL(MINUTE)
      R_DATE(8)=REAL(ISEC)
!
!-----------------------------------------------------------------------
!  Build reference-time string.
!-----------------------------------------------------------------------
!
      WRITE (TEXT,10) IYEAR, MONTH, IDAY, IHOUR, MINUTE, ISEC
 10   FORMAT (I4,'-',I2.2,'-',I2.2,1X,I2.2,':',I2.2,':',I2.2)
      REF_ATT=TEXT
      RETURN
      END FUNCTION REF_ATT

   SUBROUTINE  CONVERSION(NUMBER, YEAR, MONTH, DAY)
   IMPLICIT  NONE

   INTEGER, INTENT(IN)  :: NUMBER
   INTEGER, INTENT(OUT) :: YEAR, MONTH, DAY

   YEAR  = NUMBER / 10000
   MONTH = MOD(NUMBER, 10000) / 100
   DAY   = MOD(NUMBER, 100)
   END SUBROUTINE  CONVERSION
